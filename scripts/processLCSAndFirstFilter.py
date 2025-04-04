from .LCS import LCS
from pathlib import Path
import os

class processLCSAndFirstFilter:
    def __init__(self, outdir, outputDir_rawBlock, ratio, blockDir, sequenceDir, sourceSyntenyPath, sp, chr_shape):
        self.outdir_synteny = outdir
        self.outdir_rawBlock = outputDir_rawBlock
        self.sp = sp
        self.ratio = ratio
        self.blockDir = blockDir
        self.sequenceDir = sequenceDir
        self.sourceSyntenyPath = sourceSyntenyPath
        self.chr_shape = chr_shape

        #获取species列表
        # for i in gff_path_list:
        #     self.sp.append(i.split('/')[-1].split('.')[0])

        #获取sp之间的拷贝数比例
        # ratioDir = {}
        # with open(ratio, 'r') as rf:
        #     for line in rf:
        #         content = line.rstrip('\n')
        #         content = content.split('\t')
        #         ratioDir[content[0]] = content[1]
        # for i in self.sp:
        #     self.ratio += ratioDir[i] + ':'
        # self.ratio = self.ratio[:-1]

        if (not Path(self.outdir_rawBlock).exists()):
            os.makedirs(self.outdir_rawBlock)


    def readBlockSequence(self, filePath):
        # 读入各物种的原始block信息
        chr = []
        chr_list = []
        chr_count = 1
        fr = open(filePath, 'r')
        for line in fr:
            line = line[:-2]
            itemset = line.split(' ')[1:]
            chr.append(itemset)
            chr_list.append(str(chr_count))
            chr_count += 1
        chr_Dict = {}
        for i in range(len(chr_list)):
            chr_Dict[chr_list[i]] = chr[i]
        fr.close()
        return chr_Dict

    def readOriginSequence(self, filePath):
        # 读入原始的基因ID序列
        chr = []
        chr_list = []
        chr_count = 1
        fr = open(filePath, 'r')
        for line in fr:
            line = line[:-2]
            itemset = line.split(' ')
            chr.append(itemset)
            chr_list.append(str(chr_count))
            chr_count += 1
        chr_Dict = {}
        for i in range(len(chr_list)):
            chr_Dict[chr_list[i]] = chr[i]
        fr.close()
        return chr_Dict

    def syntenyDict(self, filePath):
        syntenyDict = {}
        with open(filePath, 'r') as sf:
            for line in sf:
                temp = line.rstrip('\n').rstrip()
                itemset = temp.split(' ')
                header = itemset[0].split(':')
                syntenyDict[header[0]] = itemset[1:]

        return syntenyDict

    def syntenyCount(self, filePath):
        syntenyDict = {}
        with open(filePath) as sf:
            for line in sf:
                temp = line.rstrip('\n').rstrip()
                itemset = temp.split(' ')
                header = itemset[0].split(':')
                if header[0] not in syntenyDict.keys():
                    syntenyDict[header[0]] = len(itemset[1:])
                else:
                    syntenyDict[header[0]] += len(itemset[1:])

        return syntenyDict
    def assambleDrimmSequence(self, blockSequence, synteny):
        # 将Block序列通过synteny文件还原为基因ID形式
        sequences = {}
        sequences_ID = {}
        blockCount = {}
        for i in blockSequence.keys():
            sequence = []
            sequence_ID = []
            for j in blockSequence[i]:
                if j.startswith('-'):
                    block = j[1:]
                    synteny_sequence = synteny[block][::-1]
                    if block not in blockCount.keys():
                        blockCount[block] = 1
                    else:
                        blockCount[block] += 1
                    for k in range(len(synteny_sequence)):
                        sequence.append(synteny_sequence[k])
                        sequence_ID.append('-' + block + '|' + str(blockCount[block]) + '|' + str(k))
                else:
                    block = j
                    synteny_sequence = synteny[block]
                    if block not in blockCount.keys():
                        blockCount[block] = 1
                    else:
                        blockCount[block] += 1
                    for k in range(len(synteny_sequence)):
                        sequence.append(synteny_sequence[k])
                        sequence_ID.append('+' + block + '|' + str(blockCount[block]) + '|' + str(k))
            sequences[i] = sequence
            sequences_ID[i] = sequence_ID
        return sequences, sequences_ID

    # 两个序列做匹配
    def matchingSequence(self, species_all_sequences, species_reassamble_sequences, species_all_sequences_name, species_reassamble_sequences_ID):
        block_range = {}
        for i in species_all_sequences.keys():
            species_reassamble_sequence = species_reassamble_sequences[i]
            # print(i)
            block_range[i] = {}
            for j in species_all_sequences[i].keys():
                if j not in species_reassamble_sequence.keys():
                    continue
                else:
                    # print(j)
                    p = LCS()
                    p.input(species_all_sequences[i][j]
                            ,species_reassamble_sequence[j])
                    direction_list, lcslength_list = p.Compute_LCS()
                    lcs = p.printOneLCS()
                    for k in lcs:
                        genename = species_all_sequences_name[i][j][k[0]]
                        ID = species_reassamble_sequences_ID[i][j][k[1]]
                        ID_split = ID.split('|')
                        block = ID_split[0][1:]
                        block_count = ID_split[1]
                        block_stand = ID_split[0][0]
                        if block not in block_range[i].keys():
                            block_range[i][block] = {}
                            block_range[i][block][block_count+'@'+block_stand] = [[j,k[0],genename]]
                        else:
                            if block_count+'@'+block_stand not in block_range[i][block].keys():
                                block_range[i][block][block_count+'@'+block_stand] = [[j,k[0],genename]]
                            else:
                                block_range[i][block][block_count+'@'+block_stand].append([j,k[0],genename])
        return block_range

    def outSynteny(self, block_range, species_all_sequences_name, species_all_sequences):
        # 输出synteny文件
        for i in block_range.keys():
            outfile = self.outdir_synteny + i + '.synteny'
            outfile_name = self.outdir_synteny + i + '.synteny.genename'
            outfile = open(outfile, 'w')
            outfile_name = open(outfile_name, 'w')
            for j in block_range[i].keys():
                block = j
                for k in block_range[i][j].keys():
                    block_count = k.split('@')[0]
                    block_stand = k.split('@')[1]
                    matching_pairs = block_range[i][j][k]
                    # print(matching_pairs)
                    matching_pairs = sorted(matching_pairs, key=lambda x: x[1])
                    chr = matching_pairs[0][0]
                    start = matching_pairs[0][1]
                    end = matching_pairs[-1][1]

                    if block_stand == '-':
                        genename = species_all_sequences_name[i][chr][start:end + 1][::-1]
                        genesequence = species_all_sequences[i][chr][start:end + 1][::-1]
                    else:
                        genename = species_all_sequences_name[i][chr][start:end + 1]
                        genesequence = species_all_sequences[i][chr][start:end + 1]
                    if len(genename) > 2:
                        outfile.write(block + ':' + str(block_count) + ':chr_' + chr + ':' + block_stand + ' ')
                        outfile_name.write(block + ':' + str(block_count) + ':chr_' + chr + ':' + block_stand + ' ' )
                        for l in genename:
                            outfile_name.write(l + ' ')
                        outfile_name.write('\n')
                        for l in genesequence:
                            outfile.write(l + ' ')
                        outfile.write('\n')
                    else: print(block + ':' + str(block_count) + ':chr_' + chr + " has " + str(len(genename)) + " anchors, drop!")
            outfile.close()
            outfile_name.close()
    
    def process_blocks(self, blocks_ratio, blockSequences, species, ratio):
        """
        Processes blocks by dynamically filtering based on species ratios and writing them into separate files.
        Args:
            blocks_ratio (dict): A dictionary mapping block IDs to their species ratios.
            blockSequences (dict): A dictionary mapping species to their block sequences.
            species (list): A list of species names.
            ratio (int): Maximum ratio value to dynamically generate target ratios (e.g., 1:1 to 1:N).
        """
        # Path for the block information file
        out_block_info_path = f"{self.outdir_synteny}/blocks.info"

        # Dynamically generate target ratios
        target_ratios = {f"1:{i}" for i in range(1, ratio + 1)}
        filters = {r: [] for r in target_ratios}

        # Write block information and filter blocks
        with open(out_block_info_path, 'w') as out_block_info:
            header = "blockID" + "\tratio(" + ":".join(species) + ")\n"
            out_block_info.write(header)

            for block_id, block_ratio in blocks_ratio.items():
                out_block_info.write(f"{block_id}\t{block_ratio}\n")
                if block_ratio in target_ratios:
                    filters[block_ratio].append(block_id)
                else:
                    print(f"{block_id} has ratio {block_ratio}, drop!!!")

        # Write filtered blocks to files
        for spec in species:
            # Prepare file paths for dynamically generated ratios
            #file_paths = {r: f"{self.outdir_rawBlock}{spec}_{r.replace(':', '_')}.block" for r in target_ratios}
            file_paths = {r: f"{self.outdir_rawBlock}{spec}_{r.split(':')[1]}.block" for r in target_ratios}
        
            # Open files and store handles in a dictionary
            output_files = {r: open(path, 'w') for r, path in file_paths.items()}
            try:
                block_sequences = blockSequences[spec]

                for chrom, sequences in block_sequences.items():
                    prefix = 's ' if self.chr_shape.lower() == 's' else 'c '
                    # Write the chromosome prefix
                    for out in output_files.values():
                        out.write(prefix)

                    # Write blocks to the appropriate files based on the filter
                    for seq in sequences:
                        block = seq[1:] if seq.startswith('-') else seq
                        for r, blocks in filters.items():
                            if block in blocks:
                                output_files[r].write(seq + ' ')
                
                    # End the line for all files
                    for out in output_files.values():
                        out.write('\n')
            finally:
                # Close all output files
                for out in output_files.values():
                    out.close()


    def excute(self):
        block_dir = self.blockDir
        species = self.sp
        #speciesRatio = self.ratio
        syntenyFile = self.sourceSyntenyPath
        sequences_dir = self.sequenceDir

        # speciesRatio = '2:2:1:1:1:2:1:1:1:1'
        # species = ['R64','castellii','gossypii','kluyveri','lactis','naganishii',
        #            'rouxii','thermotolerans','waltii','Gorden09']
        #
        # species = ['Fusarium_culmorum']
        # blockFilelist = []

        blockSequences = {}
        for i in species:
            blockSequences[i] = self.readBlockSequence(block_dir  + i + '.block')
        synteny = self.syntenyDict(syntenyFile)
        species_reassamble_sequences = {}
        species_reassamble_sequences_ID = {}
        for i in blockSequences.keys():
            # print(i)
            sequences, sequences_ID = self.assambleDrimmSequence(blockSequences[i], synteny)
            species_reassamble_sequences[i] = sequences
            species_reassamble_sequences_ID[i] = sequences_ID

        # 读取allsequence和name
        species_all_sequences = {}
        species_all_sequences_name = {}
        for i in species:
            # print(i)
            species_all_sequences[i] = self.readOriginSequence(sequences_dir + i + '.all.sequence')
            species_all_sequences_name[i] = self.readOriginSequence(sequences_dir + i + '.all.sequence.genename')

        # 两个序列做匹配
        block_range = self.matchingSequence(species_all_sequences, species_reassamble_sequences,
                                            species_all_sequences_name, species_reassamble_sequences_ID)
        self.outSynteny(block_range, species_all_sequences_name, species_all_sequences)
        # 通过synteny文件找满足条件的block，过滤不满足比例的

        blocks_ratio = {}
        block_list = []
        species_block_dict = {}
        for i in species:
            species_block_dict[i] = {}
            block_synteny_file = self.outdir_synteny + i + '.synteny'
            with open(block_synteny_file, 'r') as bsf:
                for line in bsf:
                    temp = line.rstrip('\n').rstrip()
                    itemset = temp.split(' ')
                    header = itemset[0].split(':')
                    block = header[0]
                    if block not in species_block_dict[i].keys():
                        species_block_dict[i][block] = 1
                    else:
                        species_block_dict[i][block] += 1
                    if block not in block_list:
                        block_list.append(block)


        for i in block_list:
            ratio_ = ''
            for j in species:
                if i not in species_block_dict[j].keys():
                    ratio_ += '0:'
                else:
                    ratio_ += str(species_block_dict[j][i])+':'
            ratio_ = ratio_[:-1]
            blocks_ratio[i] = ratio_

        
        out_block_info = self.outdir_synteny + '/blocks.info'
        out_block_info = open(out_block_info, 'w')
        header = "blockID" + "\tratio(" + ":".join(species) + ")\n"
        out_block_info.write(header)
        ploidy = self.ratio

        self.process_blocks(blocks_ratio, blockSequences, species, ploidy)
        ## compared with processLCSAndFirstFilter_2sp, process_blocks() replaced the final parts of that