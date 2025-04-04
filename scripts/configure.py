#!/usr/bin/env python
import os,shutil,sys
from pathlib import Path

def makeGenespaceScript(workDir, genomeIDs, ploi):
    # Ensure workDir is a valid path and has a trailing slash
    workDir = str(Path(workDir).resolve()) + "/"
    genespaceDir = workDir + "genespace"
    rScript = workDir + "run_genespace.R"

    # Ensure genomeIDs and ploi are lists
    if not isinstance(genomeIDs, list) or not isinstance(ploi, list):
        raise ValueError("genomeIDs and ploi must be lists")

    # Format lists for R script
    spec = f'c({", ".join(f'"{item}"' for item in genomeIDs)})'
    ploiStr = f'c({", ".join(f'"{item}"' for item in ploi)})'

    # R script content
    r_script_content = f'''#!/bin/env R

library(GENESPACE)

# Initialize Genespace
gpar <- init_genespace(
    wd = "{genespaceDir}",
    genomeIDs = {spec},
    ploidy = {ploiStr},
    # ignoreTheseGenomes = "ANC",  # Uncomment if needed
    minPepLen = 30,
    nCores = 8,
    maxOgPlaces = 24,
    orthofinderInBlk = FALSE
)

# Run Genespace
gpar <- run_genespace(gsParam = gpar)
'''

    # Write the content to the R script file
    with open(rScript, 'w') as file:
        file.write(r_script_content)

    return rScript  # Returning the path to the generated script

def configure(workDir,packageDir,step):
    if not workDir[-1] == "/": workDir += "/"
    if not packageDir[-1] == "/": packageDir += "/"
    genespaceDir = workDir + "genespace/"
    resultsDir = workDir + "results/"
    bedDir0 = workDir + "bed/"
    pepDir0 = workDir + "pep/"
    bedDir1 = packageDir + "bed/"
    pepDir1 = packageDir + "pep/"
    bedDir2 = genespaceDir + "bed/"
    pepDir2 = genespaceDir + "peptide/"
    meta = workDir + "meta"
    index = workDir + "index/"
    if step == "1":
        if Path(genespaceDir).exists():
            print("a floder named genespace exists, please delete it and rerun the configure ...")
        else:
            os.makedirs(genespaceDir)
            bedDir2 = genespaceDir + "bed/"
            pepDir2 = genespaceDir + "peptide/"
            os.makedirs(bedDir2)
            os.makedirs(pepDir2)
            for file_name in os.listdir(bedDir0):
                source_path = os.path.join(bedDir0, file_name)
                destination_path = os.path.join(bedDir2, file_name)
                shutil.copy2(source_path, destination_path)  # copy2 preserves metadata
            for file_name in os.listdir(pepDir0):
                source_path = os.path.join(pepDir0, file_name)
                destination_path = os.path.join(pepDir2, file_name)
                shutil.copy2(source_path, destination_path)  # copy2 preserves metadata
            for file_name in os.listdir(bedDir1):
                source_path = os.path.join(bedDir1, file_name)
                destination_path = os.path.join(bedDir2, file_name)
                shutil.copy2(source_path, destination_path)  # copy2 preserves metadata
            for file_name in os.listdir(pepDir1):
                source_path = os.path.join(pepDir1, file_name)
                destination_path = os.path.join(pepDir2, file_name)
                shutil.copy2(source_path, destination_path)  # copy2 preserves metadata
            dataPath = Path(bedDir2)
            genomeIDs = [f.stem for f in dataPath.glob("*.bed") if f.is_file()]
            ploi = ["1"] * len(genomeIDs)
            makeGenespaceScript(workDir,genomeIDs,ploi)
            print("\n====================================\nAn R script has been generated in: " + workDir + ", Please use it to run Genespace")
        folders = genespaceDir,resultsDir,meta,index,bedDir2,pepDir2
    elif step == "2":
        print("\n====================================\n\nChecking if all inputs are valid ....\n")
        # do some check here: if all the folders are here and the inputs are correct
        folders = genespaceDir,resultsDir,meta,index,bedDir2,pepDir2
        print("All inputs are ready, let's start the AGB analysis ....\n\n====================================\n")
    return folders