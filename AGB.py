#!/use/bin/env python

# get dependencies ready
import os,shutil,sys
import pandas as pd
import numpy as np
from pathlib import Path
import scripts.makeGenomeInfo as sGI
import scripts.processDrimm as spd
import scripts.processOG as spo
import scripts.modules as m
import scripts.configure as c

def main(packagedir,workdir,step):
    if not packagedir[-1] == "/": packagedir += "/"
    if not workdir[-1] == "/": workdir += "/"
    if str(step) == "1":
        info = c.configure(workdir,packagedir,step)
    elif str(step) == "2":
        info = c.configure(workdir,packagedir,step)
        drimmPath = packagedir + "software/DRIMM-Synteny"
        resultsDir = info[1]
        genome_meta = info[2]
        index_folder = info[3]
        bed_folder = info[4]
        pep_folder = info[5]
        pangenome = resultsDir + "AGB_synOG.txt"
        ## let's first clean the pangenome table from Genespace;
        clean_pangenome = m.pangenome_cleaner(pangenome)
        ## let's then parse the cleaned pangenome
        orthogroups = m.parse_pangenome(pangenome,clean_pangenome)
        # this is a tuble of the path to the outputs: *_all, *_syn, and *_synteny; *_syn will be used for the pipeline
        # Get orthogroups ready
        ortho = pd.read_csv(orthogroups[1], sep='\t', dtype = str, keep_default_na=False)
        ortho2 = ortho[ortho["interpChr"] != ""]
        
        #Step1: let's prepare inputs for the pipeline
        df_index = sGI.make_index(index_folder)
        
        DIC = {f"{row['sp']}@{row['chr']}": row['length'] for _, row in df_index.iterrows()}
        chromosome_meta = sGI.make_chromosome_meta(bed_folder,DIC)
        genomeMeta = sGI.make_genome_meta(genome_meta)
        bed_dir = resultsDir + "drimmBED/"
        if not os.path.exists(bed_dir):
            os.makedirs(bed_dir)
            sGI.make_drimmBED(bed_folder, DIC, bed_dir)
        else:
            print(bed_dir + " exists, let's generate bed file in this floder")
            sGI.make_drimmBED(bed_folder, DIC, bed_dir)

        #Step2: perform drimm-synteny to call synteny blocks
        reference = "AGB"
        queries = genomeMeta['sp'].to_list()
        genome_meta_dict = genomeMeta.set_index('sp')['chrN'].to_dict()
        chromosome_meta_dict = {lst[0]: lst[1:] for lst in chromosome_meta}
        for species in queries:
            if not species == reference:
                species_list = [reference, species]
                subdir = resultsDir + species + "/1_SynOG/" # 1. this is for two species, if we want to for multiple species, we need to change here
                info_list = m.prepare_meta(species_list, genomeMeta, chromosome_meta)
                sp,sp_ratio,sp_chr_number,long_chr_list,gff_list = info_list # this decides species pair
                # step 1, let's process cleaned pangene table to generate drimm.sequence for drimm-synteny
                if not Path(subdir).exists():
                    os.makedirs(subdir)
                    ortho3 = m.prepare_OG(ortho2,sp) # 2. here, the arguement 'sp' deciedes pair-wise or multiple-species comparasion
                    group_dir = spo.get_group_dir(ortho3,sp) # 3. and here, the arguement 'sp' deciedes pair-wise or multiple-species comparasion
                    finalGroup = spo.get_final_group(group_dir,sp_ratio)
                    spo.processSynOG(bed_dir,subdir,group_dir,finalGroup,gff_list,long_chr_list)
                else: print("SynOG for " + species + " exists, we will use this for Drimm-Synteny; If this is not correct, please delete 1_SynOG/ and rerun ...")

                # step 2, let's process drimm.sequence to generate row drimm block
                print("\n=============================================================\nlet's process orthologs to generate row blocks\n=============================================================\n")
                drimmSequence = subdir + "drimm.sequence"
                sp_ploidy = dict(zip(genomeMeta['sp'], genomeMeta['ploidy']))
                dustThreshold = sp_ploidy[species] + sp_ploidy[reference] + 1
                outPath = resultsDir + species + "/2_DrimmRaw/"
                if not Path(outPath).exists():
                    os.makedirs(outPath)
                    m.drimmSynteny(drimmSequence,drimmPath,outPath,dustThreshold)
                    print("drimm-synteny done!")
                else:
                    print("Drimm-Synteny outpath exists!!, we will check if blocks has been built ...")
                    Block_File = outPath + "blocks.txt"
                    if Path(Block_File).exists() and Path(Block_File).stat().st_size > 0:
                        print("A block file seems has been generated, we will use this for next step; If this is not the correct, please delete the block file and rerun ...")
                    else:
                        m.drimmSynteny(drimmSequence,drimmPath,outPath,dustThreshold)
                        print("drimm-synteny done!")
                
                # step 3. let's process row drimm blocks to generate final blocks
                print("\n=============================================================\nlet's process raw drimm blocks ...\n=============================================================\n")
                spd.main(info_list,species,resultsDir)
        
                print("\n=============================================================\nlet's generate cleaned blocks ...\n=============================================================\n")
                ## let's first make a combed file based on normal bed, this special bed will be used by combed_parser() and grimmBlock_parser()
                comBed = m.make_comBed(bed_folder)
                m.synBlock(resultsDir, species, comBed)
        
                print("\n=============================================================\nlet's generate data for chromosome painting ...\n=============================================================\n")
                m.chroBlock(resultsDir, species, reference, chromosome_meta_dict, genome_meta_dict)
        
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("python AGB.py packagedir workdir step")
        sys.exit(0)
    else:
        main(sys.argv[1],sys.argv[2],sys.argv[3])