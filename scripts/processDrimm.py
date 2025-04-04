# %%
import os,shutil
from pathlib import Path
from scripts.modules import prepare_meta
from scripts import processLCSAndFirstFilter as plff
from scripts import processFinalFilter as pff


def readSequence(file):
    sequence = []
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-2]
            if not line:
                break
            itemset = line.split(' ')
            sequence.append(itemset)
    return sequence

def prepare_dir(DIR):
    if DIR[-1] != "/": DIR += "/"
    drimm_split_blocks_dir = DIR + '3_DrimmBlocks/'
    raw_block_dir = DIR + '3_DrimmBlocks/tmp/'
    result_dir = DIR + '3_DrimmBlocks/finalBlocks/'
    if (not Path(drimm_split_blocks_dir).exists()):
        os.makedirs(drimm_split_blocks_dir)
    if (not Path(raw_block_dir).exists()):
        os.makedirs(raw_block_dir)
    if (not Path(result_dir).exists()):
        os.makedirs(result_dir)
    return drimm_split_blocks_dir,raw_block_dir,result_dir

def write_raw_blocks(sequence,chr_number,drimm_split_blocks_dir,sp_list):
    sp_sequences = []
    last = 0
    for i in range(len(chr_number)):
        sp_sequences.append(sequence[last:last+chr_number[i]])
        last += chr_number[i]

    for i in range(len(sp_sequences)):
        outfile = drimm_split_blocks_dir + '/' + sp_list[i] + '.block'
        outfile = open(outfile,'w')
        for j in sp_sequences[i]:
            outfile.write('s ')
            for k in j:
                outfile.write(k+' ')
            outfile.write('\n')
        outfile.close()

# to run all species-pairs in one run

#def processDrimm(info_list,working_dir,genomeMeta,chromosome_meta):
def main(info_list,name,working_dir):
    #species_list = [reference, species]
    sp_working_dir = working_dir + name + "/"
    homolog_dir = sp_working_dir + "1_SynOG/"
    if not Path(sp_working_dir).exists():
        print("no homologs found for " + name)
    else:
        block_file = sp_working_dir + "2_DrimmRaw/blocks.txt"
        sequence = readSequence(block_file)
        drimmSyntenyFile = sp_working_dir + "2_DrimmRaw/synteny.txt"
        
        #list_ = prepare_meta(species_list, genomeMeta, chromosome_meta)
        sp_list = info_list[0]
        sp_ratio = info_list[1]
        #sp_ratio = ":".join([str(e) for e in sp_ratio]) # convert sp_ratio to a string
        sp_ratio = int(sp_ratio[1])
        print(sp_ratio)
        chr_number = info_list[2]

        drimm_split_blocks_dir,raw_block_dir,result_dir = prepare_dir(sp_working_dir)
        #write_raw_blocks(sequence,chr_number)
        write_raw_blocks(sequence,chr_number,drimm_split_blocks_dir,sp_list)

        processLCSAndFirstFilter = plff.processLCSAndFirstFilter(drimm_split_blocks_dir, raw_block_dir, sp_ratio,
                                                         drimm_split_blocks_dir, homolog_dir, drimmSyntenyFile,
                                                         sp_list, 's')
        processLCSAndFirstFilter.excute()

        processFinalFilter = pff.processFinalFilter(sp_list, raw_block_dir, drimm_split_blocks_dir, result_dir, sp_ratio, 's')
        processFinalFilter.excute()
        shutil.rmtree(raw_block_dir)
