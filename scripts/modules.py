import pandas as pd
import numpy as np
import os
from itertools import combinations

def getFilterSequence(bed, group_filter_dir, long_chr_list): # add the arguement "long_chr_list" which is a list containing all long chromosomes of all species
    bed = pd.read_csv(bed,sep='\t',header=None)[[0,1]]

    chrlist = bed[0].unique().tolist()
    chrlist_long = [chr for chr in chrlist if chr in long_chr_list] # added by tao to keep only the long chromosomes

    new_sequence = []
    for i in chrlist_long: # chrlist --> chrlist_long
        split = bed.loc[bed[0] == i][1].tolist()

        chr = []
        for j in split:
            if j in group_filter_dir.keys():
                chr.append(group_filter_dir[j])
        new_sequence.append(chr)
    return new_sequence

def getAllSequence(bed, group_dir, long_chr_list): # add the arguement "long_chr_list" which is a list containing all long chromosomes of all species
    bed = pd.read_csv(bed,sep='\t',header=None)[[0,1]]
    
    chrlist = bed[0].unique().tolist()
    chrlist_long = [chr for chr in chrlist if chr in long_chr_list] # added by tao to keep only the long chromosomes

    new_sequence = []
    new_name_sequence = []
    for i in chrlist_long: # chrlist --> chrlist_long
        split = bed.loc[bed[0] == i][1].tolist()

        chr = []
        chr_name = []
        for j in split:
            if j in group_dir.keys():
                chr.append(group_dir[j])
                chr_name.append(j)

        new_sequence.append(chr)
        new_name_sequence.append(chr_name)
    return new_sequence,new_name_sequence

def outSequence(sequence,outfile):
    outfile = open(outfile,'w')
    for i in sequence:
        for j in i:
            outfile.write(str(j)+' ')
        outfile.write('\n')
    outfile.close()

def prepare_meta(species,genomeMeta,chromosome_meta):
    # let's gather information of genomes, including the species name, ploidy level and chromosome number
    filtered_df = genomeMeta[genomeMeta['sp'].isin(species)]
    sp = filtered_df['sp'].tolist()
    sp_ratio = filtered_df['ploidy'].tolist() # this decides the synteny depth
    sp_chr_number = filtered_df['chrN'].tolist()
    dic = dict(zip(filtered_df['sp'], filtered_df['chrN']))
    sp_chr = {}
    for key in dic.keys():
        if key in species:
            sp_chr[key] = dic[key]
    # print(sp_chr)
    # let's extract long chromosomes
    long_chr_list = [] # all long chromosomes in all species
    # chromosome_meta = [[sp03v5,sp03v5@a11_1_1679,sp03v5@a5_2_1622],[sp46,sp46@chr01_1_4688,sp46@chr05_2_4559]]
    for e in chromosome_meta:
        sp_name = e[0]
        #print(sp_name)
        if sp_name in sp_chr.keys():
            long_chr_number = int(sp_chr[sp_name]) # get the number of chromosomes from genome_meta
            long_chr_list_raw = e[1:long_chr_number+1] # extract the first "long_chr_number" chromosomes from chromosome_meta; this can work because chromosomes are sorted by their size (number of genes) in chromosome_meta
            for chr_name_raw in long_chr_list_raw:
                chr_name = "_".join(chr_name_raw.split("_")[0:-1]) # remove the surfix tag (_geneNumber) which was added to chromosomes by "parse_genespaceBED_for_drimmBED()"
                long_chr_list.append(chr_name)
        
    # let's prepare bed data
    gff_list = [] #should be bed in format
    for i in sp:
        gff_list.append(i + '.bed')
    
    return sp,sp_ratio,sp_chr_number,long_chr_list,gff_list

def prepare_OG(ortho,species,OG_type='synOG'):
    #ortho = orthogroups
    #ortho = pd.read_csv(orthogroups, sep='\t', dtype = str, na_values='NA')
    if OG_type == "synOG": #works with synOG
        ortho = ortho.apply(lambda x: x.astype(str))
        ortho = ortho.apply(lambda x: x.str.replace('|', ','))
        ortho = ortho.apply(lambda x: x.str.split(',').apply(lambda y: ', '.join([item.split('@')[1] if '@' in item else item for item in y]))) # this manipulate the synOG table to make it same with orthofinder OG table
        #print(ortho)
        list_ = ["pgID"] + species # narrow down the df to target genomes
    else: # work with orthofinder derived OGs
        list_ = ["HOG"] + species # narrow down the df to target genomes

    ortho = ortho[list_]
    #ortho = ortho.drop(ortho.columns[[1, 2]], axis=1) # drop the 2nd and 3rd columns
    ortho = ortho.fillna('')
    #columns = ortho.columns.tolist() # what does this do?
    ortho = np.asarray(ortho)
    return ortho

def standarlize_drimmblocks(genome_list, drimmBlocks, grimmBlocks):
    # for a list of species/genomes (genome_list), convert the drimm blocks (drimmBlocks) to grimm blocks (grimmBlocks)for further analysis by grimm
    # let's gather all drimm genomes (namely final.block files from process_drimm), including extant ones and their ancestors, put them in a dictionary
    dic_ = {}
    for genome in drimmBlocks:
        genome_name = os.path.basename(genome).split(".")[0]
        chromosomes = []
        for line in open(genome,"r"):
            if line.startswith("s"): # in case there are blank lines
                line = line.strip().split(" ")
                line = "~".join(line[1:])
                chromosomes.append(line)
        dic_[genome_name] = chromosomes # {genome1:[g1~g2,g3~g4,g5~g6],genome2:[]}
    
    # let's use the dictionary generated above to index the genes, 1 to x
    dic__ = {}
    i = 1
    for genomes in dic_.keys():
        chromosomes = dic_[genomes]
        for chromosome in chromosomes:
            genes = chromosome.split("~")
            for gene in genes:
                gene = gene.strip("-")
                if not gene in dic__.keys():
                    dic__[gene] = str(i)
                    i += 1
    
    # let's make use of the index and reformat the genomes
    with open(grimmBlocks, "w") as f:
        for genome in dic_.keys():
            print(genome)
            if genome in genome_list: # a list of the genomes we want to include, sometimes we want only analysis some geneomes, for example drimm consider ancestor to descent, mgra takes only extant genomes
                print(genome + " is needed, let's process it....")
                f.write(">" + genome + "\n")
                chromosomes = dic_[genome]
                for chromosome in chromosomes:
                    genes = chromosome.split("~")
                    genes_new = []
                    for gene in genes:
                        if gene.startswith("-"):
                            gene = gene.strip("-")
                            gene_new = "-" + dic__[gene]
                        else:
                            gene_new = dic__[gene]
                        genes_new.append(gene_new)
                    chromosome_new = " ".join(genes_new) + " $"
                    f.write(chromosome_new + "\n")
                f.write("\n") # write a blank line between genomes
    f.close()

def processGenenumber(sp, resultDir, ratio):
    block_len = {}
    for i in sp:
        with open(resultDir + '/' + i + '_' + str(ratio) + '.final.synteny', 'r') as sf:
            for line in sf:
                temp = line
                temp = temp.rstrip('\n').rstrip()
                block = temp.split(' ')[0].split(':')[0]
                length = len(temp.split(' ')[1:])
                # print(temp.split(' ')[1:])

                if block not in block_len.keys():
                    block_len[block] = length
                else:
                    if block_len[block] < length:
                        block_len[block] = length

    with open(resultDir + '/blockindex.genenumber_' + str(ratio),'w') as f:
        f.write('blockID\tblockLength\n')
        for i,j in block_len.items():
            f.write(i + '\t' + str(j) + '\n')
    return

def block_parser(block_file, species, dic):
    list__ = []
    for line in open(block_file, "r"):
        line = line.strip()
        block_id = line.split(" ")[0] # "4099:1:chr_1:-"
        block_member = line.split(" ")[1:]
        #block_left_member = line.split(" ")[1]
        #block_right_member = line.split(" ")[-1]
        coo = []
        ord = []
        for member in block_member:
            member = species + "@" + member
            if not member in dic.keys():
                print(member + " not found in combed")
            else:
                member_new = dic[member] # [chr, start, end, nGeneOnChr, ord]
                coo.append(int(member_new[1]))
                coo.append(int(member_new[2]))
                ord.append(int(member_new[4]))
                chrLength = member_new[3]
        coo.sort()
        boundary_coo = str(coo[0]) + "-" + str(coo[-1])
        block_length_bp = coo[-1] - coo[0]
        block_length_gene = len(block_member)
        ord.sort()
        boundary_ord = str(ord[0]) + "-" + str(ord[-1])
        ord = [str(e) for e in ord]
        member_ord = ','.join(ord)
        list_ = [block_id, block_length_gene, boundary_ord, boundary_coo, block_length_bp, chrLength, member_ord]
        list__.append(list_) # make a list of list
    columns = ["block_id", "block_length_gene", "boundary_geneOrd", "boundary_coo", "block_length_bp", "chrLength", "member_ord"]
    block = pd.DataFrame(list__, columns=columns) # make a dataframe from the list of list
    return block

# let's first make a combed file based on normal bed, this special bed will be used by combed_parser() and grimmBlock_parser();
def make_comBed(beddir):
    bed_data_frames = []
    bedFiles = [os.path.join(beddir, filename) for filename in os.listdir(beddir) if filename.endswith('.bed')]
    for bed_file in bedFiles:
        print("processing " + os.path.basename(bed_file))
        df = pd.read_csv(bed_file, sep='\t', names=["chr", "start", "end", "id"])
        df['genome'] = os.path.splitext(os.path.basename(bed_file))[0]
        bed_data_frames.append(df)
    bed = pd.concat(bed_data_frames, ignore_index=True)
    #bed['chrName'] = bed['chr'].str.extract('(\d+)').astype(float) 
    # group the records according to 'genome' and 'chr' <groupby(['genome', 'chr'])>, then count the records in each 'genome' group <['genome'].transform('count')>, asign the count to new column 'nGeneOnChrom'
    bed['nGeneOnChr'] = bed.groupby(['genome', 'chr'])['genome'].transform('count')
    # sort the records by 'genome', 'chr', 'nGeneOnChrom', 'start' in either ascending or descending order
    bed = bed.sort_values(by=['genome', 'chr', 'nGeneOnChr', 'start'], ascending=[True, True, False, True])
    # add a column with data indicate the order of the gene in the chromsome
    bed['ord'] = bed.groupby(['genome', 'chr']).cumcount() + 1
    bed = bed.reset_index(drop=True)
    #bed['ord'] = range(1, len(bed) + 1)
    return bed

# let's read in the combed file and parse it
def combed_parser(combed):
    '''
    #if combed is a txt file, do the following:
    dict = []
    for line in open(combed, "r"):
        line = line.strip()
        line = line.split("\t")
        species_gene = line[4] + "@" + line[3] # 'sp21@LSAT1v11_C40016059-RA'
        feature = [line[0], line[1], line[2], line[5], line[6]] # [chr, start, end, nGeneOnChr, ord]
        dict[species_gene] = feature
    '''
    #if combed is a dataframe, do the following:
    dict = {f"{row['genome']}@{row['id']}": [row['chr'], row['start'], row['end'], row['nGeneOnChr'], row['ord']]
               for _, row in combed.iterrows()}
    return dict

def grimmBlock_parser(blockFiles,dic,output_path):
    for block_file in blockFiles:
    #species = os.path.splitext(os.path.basename(block_file))[0] # this split the file into name and extension (split in last dot)
        species = os.path.basename(block_file).split(".")[0]
        if species[-2] == "_":
            print("processing " + species + " blocks ...")
            species_ = species.split("_")[0]
            block_pd = block_parser(block_file, species_, dic)
            print("writing " + species + " blocks ...")
            block_pd.to_csv(output_path + species + '_synBlock.txt', sep='\t', header=True, index=False)
            print(species + " done")
        else:
            print("processing " + species + "...")
            block_pd = block_parser(block_file, species, dic)
            print("writing " + species + " blocks ...")
            block_pd.to_csv(output_path + species + '_synBlock.txt', sep='\t', header=True, index=False)
            print(species + " done")

def get_block_info(reference,target):
    '''
    1. let's get the information of origin of blocks from ancestor (or reference), namely which reference chromosome is the block from ?
    2. reference is the final.synteny.genename from sp32_synBlocks.txt (output from process_drimm_blocks.py)
    '''
    dic1 = {}
    dic2 = {}
    for line in open(reference, "r"):
        if not line.startswith("block_id"):
            blockID = line.split(":")[0]
            chrID = line.split(":")[2]
            dic1[blockID] = chrID # 1576:chr_1

    for line in open(target, "r"):
        if not line.startswith("block_id"):
            line = line.strip()
            line = line.split("\t")
            blockHeader = line[0] # 1576:1:chr_1:+
            #print(blockHeader)
            blockID = blockHeader.split(":")[0] # 1576
            blockSubgenome = blockHeader.split(":")[1] # 1
            blockChromosome = blockHeader.split(":")[2] # chr_1
            blockBoundary_gene = line[2]
            blockBoundary_bp = line[3]

            dic2[blockID + "_" + blockSubgenome + "_" + blockChromosome] = [blockBoundary_gene,blockBoundary_bp]
            # 1576_1_chr_1:,1-10,1-1000]
            # 1576_2_chr_1:[50-100,10000-20000]
    return dic1,dic2

def get_chromosome_info(target_chr,dic1,dic2):
    c = 0
    coo_gene = []
    coo_bp = []
    for line in open(target_chr, "r"):
        line = line.strip()
        c += 1
        chromosome = "chr_" + str(c)
        line = line.split(" ")[1:]
        list_ = []
        for block in line: # 1576
            block = block.strip("-") # but if we want to keep the information of block direction, how do we do?
            for i in ["_1","_2","_3"]:
                blockID = block + i + "_" + chromosome # 1576_1_chr_1
                #print(blockID)
                if blockID in dic2.keys():
                    list_.append(blockID)
        #print(list_)
        
        list__ = list(set(list_)) # we remove redundant elements
        #print("\t".join(list__))
        list___ = {}
        for e in list__:
            list___[e] = int(dic2[e][0].split("-")[0])
        sorted_dic = sorted(list___, key=list___.get) # sort the keys (blockID) based on their corresponding values (left boundary)
        sorted_blockID = list(sorted_dic) # we get a new list of blocks, which is same with the original chromosomes but with blocks distinguished in _1, _2 and _3
        print(sorted_blockID)

        # now let's get 
        for e in sorted_blockID:
            block = e.split("_")[0]
            colour = dic1[block] # chr_1
            coo_gene_left = dic2[e][0].split("-")[0]
            coo_gene_right = dic2[e][0].split("-")[1]
            coo_bp_left = dic2[e][1].split("-")[0]
            coo_bp_right = dic2[e][1].split("-")[1]
            #coo_gene.append([chromosome,e,coo_gene_left,coo_gene_right,colour])
            #coo_bp.append([chromosome,e,coo_bp_left,coo_bp_right,colour])
            coo_gene.append([chromosome,"0",coo_gene_left,coo_gene_right,colour])
            coo_bp.append([chromosome,"0",coo_bp_left,coo_bp_right,colour])
    return coo_gene,coo_bp

def drimmSynteny(drimmSequence,drimmpath,outPath,dustThreshold,cycleLength='20'):
    drimmout = outPath.replace("&", r"\&")
    drimmInput = drimmSequence.replace("&", r"\&")
    run_drimmSynteny = 'mono ' + drimmpath + ' ' + drimmInput + ' ' + drimmout + ' ' + cycleLength + ' ' + str(dustThreshold)
    print(run_drimmSynteny)
    os.system(run_drimmSynteny)

def synBlock(workdir, species, comBed):
    blockDir = workdir + species + "/3_DrimmBlocks/finalBlocks/"
    output_path = blockDir + "../../5_SynBlocks/"
    if not os.path.exists(output_path): os.mkdir(output_path)
    dic = combed_parser(comBed)
    blockFiles = [os.path.join(blockDir, filename) for filename in os.listdir(blockDir) if filename.endswith('.synteny.genename')]
    grimmBlock_parser(blockFiles,dic,output_path)

def chroBlock(workdir, species, reference, chromosome_meta_dict, genome_meta_dict):
    blockDir = workdir + species + "/5_SynBlocks/"
    blocks = [file for file in os.listdir(blockDir) if file.endswith("synBlock.txt")]
    for block in blocks:
        target_block_name = block.split("_")[0]
        target_block_type = block.split("_")[1]
        reference_block = blockDir + reference + "_" + target_block_type + "_synBlock.txt"
        target_block = os.path.join(blockDir, block)
        target_chr = blockDir + "../3_DrimmBlocks/finalBlocks/" + target_block_name + "_" + target_block_type +".final.block"
        
        dic1,dic2 = get_block_info(reference_block,target_block)
        coo_gene,coo_bp = get_chromosome_info(target_chr,dic1,dic2)

        columns = ["chr","haplotype", "start", "end", "ancestral_group"]
        output_path = blockDir + "../6_chromPainting/"
        if not os.path.exists(output_path): os.mkdir(output_path)

        r1 = pd.DataFrame(coo_gene, columns=columns) # make a dataframe from the list of list
        r2 = pd.DataFrame(coo_bp, columns=columns) # make a dataframe from the list of list
        out_put1 = output_path + target_block_name + "_" + target_block_type +"_cooGene.txt"
        out_put2 = output_path + target_block_name + "_" + target_block_type +"_cooBp.txt"
        r1.to_csv(out_put1, sep='\t', header=True, index=False)
        r2.to_csv(out_put2, sep='\t', header=True, index=False)

        chromosome_list = chromosome_meta_dict[species]
        chr_info_list = []
        i = 1
        for chromosome in chromosome_list:
            if i <= int(genome_meta_dict[species]): # make sure only take the long chromosomes as specified in genome_meta
                drimm_chr = "chr_" + chromosome.split("_")[-2]
                genome_chr = "_".join(chromosome.split("_")[0:-2])
                length = chromosome.split("_")[-1]
                chr_info = [drimm_chr,genome_chr,"0",length]
                chr_info_list.append(chr_info)
                i += 1
        df = pd.DataFrame(chr_info_list, columns=['drimm_chr', 'genome_chr', 'start','end'])
        df = df.sort_values(by='genome_chr')
        df.to_csv(output_path+species+".txt",sep="\t",index=False)

def prefix_cell(cell, header):
    if pd.notna(cell) and cell.strip():  # Check for non-empty and non-whitespace cells
        items = cell.split("|")
        items = [header + item for item in items]
        return "|".join(items)
    else:
        return "NA"

def pangenome_cleaner(pangenome):
    table = []
    f1 = open(pangenome, "r")
    for line in f1:
        # read text data into a data table
        line = line.strip()
        line = line.split("\t")
        e = [line[0],line[1],line[4]] + line[9:]
        table.append(e)
    
    df = pd.DataFrame(table)
    species_code = df.iloc[0] + "@" # extract species code from the first raw, and prefix it to genes
    df.iloc[1:, 3:] = df.iloc[1:, 3:].apply(lambda col: col.apply(lambda cell: prefix_cell(cell, species_code[col.name])))
    return df

# 2. read in the modified pangenome.txt data and converst it to nonSynteny,tandem and synteny datatable, and also datatable for SynNet package
# this script is different from 2_pangenome_parser.py by categorize genes into syntenic genes and all (syntenic, non-syntenic and tandem)
def parse_pangenome(pangenome,clean_pangenome):
    prefix = os.path.splitext(os.path.basename(pangenome))[0]
    print(prefix)
    dir_path = os.path.dirname(pangenome)
    all_ = dir_path + "/" + prefix + "_all.txt"
    syn = dir_path + "/" + prefix + "_syn.txt"
    synnet = dir_path + "/" + prefix + "_synnet.txt"
    with open(all_,"w") as f2, open(syn,"w") as f3, open(synnet,"w") as f4:
        list_ = []
        for index, row in clean_pangenome.iterrows():
            res = []
            res_new = []
            list_all = []
            list_syn = []
            list_synnet = []
            character = "+*"
            pgid = str(row.iloc[0]) # extract the pangene ID which is in the first coloum, if pangenome eas read in by pd.read_csv(), we need to convert int to str
            #note: if you don't specify the delimiter in string.split(), the split() will take mutile "\t" as a single "\t"
            list_line = row.iloc[1:].tolist() # extract the member of the given pangene from all species, and store them in list;
            #print(list_line)
            for element in list_line:
                if element == "NA": # no gene
                    list_all.append(element)
                    #list_tandem.append(element)
                    list_syn.append(element)
                elif "|" in element: # multiple genes (>=2)
                    member = element.split("|")
                    list1 = []
                    list2 = []
                    for e in member:
                        if not e.endswith("*") and not e.endswith("+"):
                            list1.append(e) # syntenic members, the associated tandems(+) are excluded
                        else:
                            e = e.strip(character)
                            list2.append(e) # other members, including syntenic-associated tandems(+) and non-syntenic members
                    if list1 == []: # all genes are non-syntenic (either * or +), means list2 is not empty
                        list1.append("NA")
                    else: # list1 is not empty, means that there are at least one syntenic genes
                        list2 = list1 + list2 # add syntenic genes to non-synetnic list
                    list_syn.append("|".join(list1))
                    list_all.append("|".join(list2))
                elif element.endswith("*"): # only one non-syntenic gene
                    element = element.strip(character)
                    list_all.append(element)
                    list_syn.append("NA")
#                elif element.endswith("+"): # this is not needed, because if there is "+", "|" must exist
#                    list_tandem.append(element)
#                    list_nonSyn.append("NA")
#                    list_syn.append("NA")
                else: # only one syntenic gene
                    list_syn.append(element)
                    list_all.append(element)

            # write results
            if pgid == "pgID": # write the header line
                f2.write(pgid + "\t" + "\t".join(list_line) + "\n")
            else:
                f2.write(pgid + "\t" + "\t".join(list_all) + "\n")
            f3.write(pgid + "\t" + "\t".join(list_syn) + "\n") # the header line has been included, see line 32     
                                    
            # generate non-redundant syntenic gene pairs of all genes stored in list_syn per chromosome
            chromosome = list_line[0]
            if not chromosome in list_ and not pgid == "pgID" and not chromosome == "":
                list_.append(chromosome)
                line = "\t".join(list_line[2:])
                line = line.replace("|", "\t") # split genes from same species as independent elements
                gene_list = line.split("\t")
                for element in gene_list:
                    if not element.endswith("*") and not element.endswith("+") and element != "NA": # this ensures that non-syntenic genes and syntenic tandem duplicates were removed
                        list_synnet.append(element)
                pairs = list(combinations(list_synnet, 2))
                for element in pairs:
                    string = "\t".join(element)
                    res.append(string)

                # remove gene pairs from same species
                for e in res:
                    gene1 = e.split("\t")[0]
                    #id1 = gene1.split("_")[0:3]
                    name1 = gene1.split("_")[0]
                    gene2 = e.split("\t")[1]
                    #id2 = gene2.split("_")[0:3]
                    name2 = gene2.split("_")[0]
                    if name1 != name2: # this will exclude pairs from the same species, in my case sp21 has 3 varieties and all will be removed if a pair has any of them
                        res_new.append(e)
                        #res.remove(e) # there is probablem with this, the iteration will stop after first match
                if len(res_new) >= 1:
                    res_new = [chromosome + "\t" + element for element in res_new]
                    #f4.write("### the clusters on " + chromosome + "\n")
                    f4.write("\n".join(res_new) + "\n")
            elif chromosome in list_ and not pgid == "pgID" and not chromosome == "":
                #print(list_line)
                line = "\t".join(list_line[2:])
                line = line.replace("|", "\t") # split genes from same species as independent elements
                gene_list = line.split("\t")
                for element in gene_list:
                    if not element.endswith("*") and not element.endswith("+") and element != "NA": # this ensures that non-syntenic genes and syntenic tandem duplicates were removed
                        list_synnet.append(element)
                pairs = list(combinations(list_synnet, 2))
                for element in pairs:
                    string = "\t".join(element)
                    res.append(string)

                # remove gene pairs from same species
                for e in res:
                    gene1 = e.split("\t")[0]
                    #id1 = gene1.split("_")[0:3]
                    name1 = gene1.split("_")[0]
                    gene2 = e.split("\t")[1]
                    #id2 = gene2.split("_")[0:3]
                    name2 = gene2.split("_")[0]
                    if name1 != name2: # this will exclude pairs from the same species, in my case sp21 has 3 varieties and all will be removed if a pair has any of them
                        res_new.append(e)
                        #res.remove(e) # there is probablem with this, the iteration will stop after first match
                if len(res_new) >= 1:
                    res_new = [chromosome + "\t" + element for element in res_new]
                    f4.write("\n".join(res_new) + "\n")
    f2.close()
    f3.close()
    f4.close()
    return all_,syn,synnet