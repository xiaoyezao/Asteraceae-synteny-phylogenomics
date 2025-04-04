#%%
import os
import pandas as pd

'''
parse_genespaceBED_for_drimmBED(): as name says, manipulate bed files to fit into drimm-synteny:
    1) prefix genome name to chromosome name, --> "genome@chromosome"
    2) count gene number on chromosomes and sort by chromosome size (gene number), chr_name then start; this makes sure that longer chromosomes come eariler in order
    3) surfix a uniq number to each chromosome, from 1 to "total chromosome number"
    4) gather the information of all chromosomes, [chr_order_size, ]
    5) reformat the bed and output it
'''
def make_index(index_folder):
    ref_index = [["AGB","a01","20946919","26","60"],
                 ["AGB","a02","15478902","21296087","60"],
                 ["AGB","a03","21342526","37032997","60"],
                 ["AGB","a04","56741044","58731258","60"],
                 ["AGB","a05","64352198","116418013","60"],
                 ["AGB","a06","20187645","181842774","60"],
                 ["AGB","a07","24478091","202366906","60"],
                 ["AGB","a08","31761148","227252992","60"],
                 ["AGB","a09","33950554","259543519","60"],
                 ["AGB","a10","18939360","294059943","60"],
                 ["AGB","a11","69615488","313314986","60"],
                 ["AGB","a12","39018945","384090760","60"],
                 ["AGB","a13","36774632","423760048","60"],
                 ["AGB","a14","28095866","461147618","60"],
                 ["AGB","a15","16859767","489711776","60"],
                 ["AGB","b01","13279672","506852566","60"],
                 ["AGB","b02","16809921","520353592","60"],
                 ["AGB","b03","16476605","537443705","60"],
                 ["AGB","b04","39003322","554194947","60"],
                 ["AGB","b05","22078130","593848351","60"],
                 ["AGB","b06","52213507","616294476","60"],
                 ["AGB","b07","17192078","669378235","60"],
                 ["AGB","b08","42167297","686856874","60"],
                 ["AGB","b09","16726989","729726986","60"],
                 ["AGB","b10","32304678","746732786","60"],
                 ["AGB","b11","61467173","779575903","60"],
                 ["AGB","b12","42130764","842067556","60"],
                 ["AGB","b13","51440046","884900527","60"],
                 ["AGB","b14","23193114","937197935","60"],
                 ["AGB","b15","6922114","960777628","60"],
                 ["AGB","c01","8861220","967815137","60"],
                 ["AGB","c02","25588691","976824070","60"],
                 ["AGB","c03","27847782","1002839266","60"],
                 ["AGB","c04","29935135","1031151204","60"],
                 ["AGB","c05","26392843","1061585284","60"],
                 ["AGB","c06","44775773","1088418034","60"],
                 ["AGB","c07","33775477","1133940096","60"],
                 ["AGB","c08","45681754","1168278524","60"],
                 ["AGB","c09","23619799","1214721667","60"],
                 ["AGB","c10","26172618","1238735157","60"],
                 ["AGB","c11","46279443","1265344013","60"],
                 ["AGB","c12","56011074","1312394808","60"],
                 ["AGB","c13","39460313","1369339427","60"],
                 ["AGB","c14","61053713","1409457439","60"],
                 ["AGB","c15","29266455","1471528741","60"]]
    sliced_data = [lst[:3] for lst in ref_index]
    df_ref = pd.DataFrame(sliced_data, columns=['sp', 'chr', 'length'])
    df_list = [df_ref]
    for index in os.listdir(index_folder):
        if index.endswith(".fai"):
            sp = os.path.splitext(index)[0]
            index_file = index_folder + index
            df_target = pd.read_csv(index_file, sep='\t', names=["chr", "length", "na1", "na2","na3"])
            df_target = df_target[["chr", "length"]]
            df_target["sp"] = sp
            df_list.append(df_target)
    result_df = pd.concat(df_list, ignore_index=True)
    return result_df

'''
def make_genome_meta(sp,ploid,chrN):
    list1 = [["AGB",1,45]]
    list2 = [sp,int(ploid),int(chrN)] # genome/species name, ploid status vs. AGB, chromosome number
    list1.append(list2)
    df = pd.DataFrame(list1,columns=['sp', 'ploid', 'chrN'])
    return df
'''

def make_genome_meta(meta):
    list1 = [["AGB",1,45]]
    list2 = []
    for line in open(meta, "r"):
        line = line.strip()
        if not line.startswith('#'):
            line = line.split()
            list2.append([line[0], int(line[1]), int(line[2])])
    list = list1 + list2
    df = pd.DataFrame(list,columns=['sp', 'ploidy', 'chrN'])
    return df

def make_chromosome_meta(BED_folder,DIC):
    chromosome_meta_dict = []
    for file in os.listdir(BED_folder):
        if file.endswith(".bed"):
            bed = BED_folder + file
            genome_name = os.path.splitext(file)[0]
            chr_info,df = parse_genespaceBED_for_drimmBED(bed,DIC)
            chr_size = ["_".join(c.split("_")[0:-1]) for c in chr_info] # remove the chr length (e.g., "sp03v1@a11_1_1679_69615488")
            list_ = [genome_name] + chr_size
            chromosome_meta_dict.append(list_)
    return chromosome_meta_dict

def make_drimmBED(BED_folder, DIC, bed_dir):
    if not BED_folder[-1] == "/": BED_folder += "/"
    for file in os.listdir(BED_folder):
        if file.endswith(".bed"):
            genome_name = os.path.splitext(file)[0]
            bed = BED_folder + file
            new_bed = bed_dir + genome_name + ".bed"
            if not os.path.exists(new_bed):
                print(new_bed + " is new, adding to database ...")
                chr_info,df = parse_genespaceBED_for_drimmBED(bed,DIC)
                df.to_csv(new_bed, header=None, index=None, sep='\t')
            else:
                print(new_bed + " exists!")

def get_chromosome_lenght(index_data):
    # gather information of chromosome_lenght from genome index(.fai)
    dic_ = {}
    for line in open(index_data, "r"):
        species = line.split("\t")[0]
        chromosome = line.split("\t")[1]
        length = line.split("\t")[2]
        dic_[species + "@" + chromosome] = length
    return dic_

def parse_genespaceBED_for_drimmBED(BED,DIC):
    df = pd.read_csv(BED, sep='\t', names=["chr", "start", "end", "gene"])
    genome_name = os.path.splitext(os.path.basename(BED))[0]
    df["chr"] = df["chr"].map(lambda x: genome_name + "@" + x)
    #df.sort_values(by=['chr', 'start'], inplace=True)

    gene_counts = df['chr'].value_counts().reset_index() # calculate number of genes on chromosome
    gene_counts.columns = ['chr', 'gene_count'] # make new columns "chr gene_number"
    df = pd.merge(df, gene_counts, on='chr') # merge df and gene_counts use "chr"
    df['length_bp'] = df['chr'].map(DIC)
    df.sort_values(by=['gene_count', 'chr', 'start'], ascending=[False, True, True], inplace=True) # be careful with df.sort_values(), by default the sorting is performed lexicographically, means that "10" would come before "2" when the column is string or mixed

    df["chr_suf"] = df.groupby("chr", sort = False).ngroup() + 1 # asign a number to the same chr in the order that is in the original df (sort = False)
    df["chr_suf"] = df["chr_suf"].astype(str)
    df["chr"] = df["chr"] + '_' + df["chr_suf"] # surfix the number to the chr name

    df["chr_size"] = df["chr"] + "_" + df["gene_count"].astype(str) + "_" + df['length_bp'].astype(str) # create a new column with "chr_geneCount_chrLength" as value
    chr_size = df["chr_size"].unique().tolist() # put chr information to a list

    df = df[['chr', 'gene', 'start', 'end']]
    #df.to_csv(genome_name + '.bed', header=None, index=None, sep='\t')

    return chr_size,df

# %%
