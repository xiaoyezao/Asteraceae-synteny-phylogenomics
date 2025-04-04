#%%
import os
from pathlib import Path
import pandas as pd
import numpy as np
#from modules_paper import prepare_OG,prepare_meta,getFilterSequence,getAllSequence,outSequence
import scripts.modules as m

def get_group_dir(ortho,sp):
    group_dir = {}
    for i in ortho:
        group = i[0]
        group_dir[group] = {}
        species = i[1:]
        for j in range(len(species)):
            genes = species[j].split(', ')
            if genes[0] == '':
                group_dir[group][sp[j]] = []
            else:
                group_dir[group][sp[j]] = genes
    return group_dir

def get_final_group(group_dir,sp_ratio):
    rate_dir = {}
    finalGroup = {}
    for i in group_dir.keys():
        rate_list = []
        for j in group_dir[i].keys():
            rate_list.append(len(group_dir[i][j]))
        ok = 1
        for j in range(len(rate_list)):
            if rate_list[j] > sp_ratio[j] or rate_list[j] == 0:
                ok = 0
        if ok == 0:
            continue
        else:
            rate = ''
            for j in rate_list:
                rate += str(j) + ':'
            rate = rate[:-1]
            finalGroup[i] = group_dir[i]
            if rate not in rate_dir.keys():
                rate_dir[rate] = 1
            else:
                rate_dir[rate] += 1
    print('gene rate')
    for i in rate_dir.keys():
        print(i + '\t' + str(rate_dir[i]))
    return finalGroup

def processSynOG(bed_dir,outdir,group_dir,finalGroup,gff_list,long_chr_list):
    if not outdir[-1] == "/": outdir += "/"
    outfile = outdir + 'group.xls'
    count = 1
    outfile = open(outfile,'w')
    outfile.write('gene\tgroup\n')

    outfile_filter = outdir + 'filter_group.xls'
    outfile_filter = open(outfile_filter,'w')
    outfile_filter.write('gene\tgroup\n')

    for i in group_dir.keys():
        for j in group_dir[i].keys():
            for k in group_dir[i][j]:
                outfile.write(k+'\t'+str(count)+'\n')
        if i in finalGroup.keys():
            for j in finalGroup[i].keys():
                for k in finalGroup[i][j]:
                    outfile_filter.write(k + '\t' + str(count) + '\n')
        count += 1
    outfile.close()
    outfile_filter.close()

    group = pd.read_csv(outdir +'group.xls',sep='\t')
    group = np.asarray(group)
    group_dir = {}
    for i in group:
        group_dir[i[0]] = i[1]

    group_filter = pd.read_csv(outdir + 'filter_group.xls',sep='\t')
    group_filter = np.asarray(group_filter)
    group_filter_dir = {}
    for i in group_filter:
        group_filter_dir[i[0]] = i[1]

    sample_sequence_files = outdir + 'drimm.sequence'
    sample_sequence_files = open(sample_sequence_files, 'w')

    for i in gff_list:
        print("let's write syntenic data for " + i.split(".")[0] + "..." )
        gff = bed_dir + i
        sequence,sequence_name = m.getAllSequence(gff, group_dir, long_chr_list) # add long_chr_list
        filter_sequence = m.getFilterSequence(gff, group_filter_dir, long_chr_list) # add long_chr_list
        item = i.split('.')
        outfile = outdir + item[0] +'.sequence'
        m.outSequence(filter_sequence, outfile)
        outallfile = outdir + item[0] +'.all.sequence'
        m.outSequence(sequence, outallfile)
        outallfilename = outdir + item[0] + '.all.sequence.genename'
        m.outSequence(sequence_name, outallfilename)

        for j in filter_sequence:
            for k in j:
                sample_sequence_files.write(str(k) + ' ')
            sample_sequence_files.write('\n')

    sample_sequence_files.close()