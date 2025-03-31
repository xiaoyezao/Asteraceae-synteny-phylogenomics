# Asteraceae-synteny-phylogenomics
<img width="743" alt="Screenshot 2025-01-14 at 21 41 59" src="https://github.com/user-attachments/assets/b35bc922-0a56-41d3-b7f3-d93d1c773560" />

## Introduction
A pipeline for Asteraceae synteny-phylogenomic analysis
>
The objective of this pipeline is to apply the synteny-phylogenomic framework built in our [paper](https://www.biorxiv.org/content/10.1101/2025.01.08.631874v1) to study Asteraceae genome evolution in context of ancient genome triplication.
>
1. we use Genespace or Mcscan to map genomes onto the 15*3 Asteraceae Genome Blocks (AGB) that generated based on the comparasion of an diploid outgroup and the paleohexaploid ingroup, from which we obtained syntenic genes as anchors
2. we use the anchors and drimm-synteny to scan the genomes to identify conserved synteny segments
3. with the synteny segments recovered, we can represent any Asteraceae genomes using the 15*3 AGB system
4. we can characterize the genome rearrangements from AGBs to focal Asteraceae genomes, and any other genome pairs as well
5. we can quantify the gene fractionation in subgenomes 

## Installation and Dependencies
### The pipeline is designed in Python3, and organised in [Jupyter Notebook](https://jupyter.org/). We have only tested it on Mac and Linux, but if all dependencies are available it should also work on Windows.
```sh
# download the package to your desired folder
git clone https://github.com/xiaoyezao/Asteraceae-synteny-phylogenomics
cd Asteraceae-synteny-phylogenomics
```
### In addition, several other softwares are needed

[Genespace](https://github.com/xiaoyezao/GENESPACE)
```sh
# I recommend making a new environment and install all softwares and library dependencies in the same Conda environment.
conda create -n AGB
conda activate AGB

conda install -c bioconda orthofinder=2.5.5
conda install bioconda::mcscanx

conda install r-data.table r-dbscan r-R.utils r-devtools
conda install bioconductor-Biostrings bioconductor-rtracklayer

devtools::install_github("xiaoyezao/GENESPACE", upgrade = F)
# if any R dependencies are still missing, try to install from Conda
```
>
DRIMM-Synteny (a pre-built executable is available in the software folder)
   ```sh
   # install mono if you don't have it on your computer
   conda install mono -c conda-forge
   # test DRIMM-Synteny
   mono DRIMM-Synteny.exe [arguements]
   
   #if this doesn't work, compile your own DRIMM-Synteny from DRIMM-Synteny.cs (available in the software folder)
   mcs DRIMM-Synteny.cs -out:DRIMM-Synteny
   # test
   mono DRIMM-Synteny --help
   ```

## Usage
#### 1. Input
Please prepare the following input data:
1) bed folder which contains all bed files
2) pep folder which contains all proteome files
3) index folder which contains all genome index files
```sh
/workingDirectory
└─ peptide
    └─ genome1.fa
    └─ genome2.fa
└─ bed
    └─ genome1.bed
    └─ genome2.bed
└─ index
    └─ genome1.index
    └─ genome2.index

```
#### 2. Call homologous groups
To call homology, Orthofinder, Mcscan and GENESPACE can be used, however, currently we found that GENESPACE works best, here we will do a demo using GENESPACE.
>
```R

wd <- "~/path/to/genespace/workingDirectory"
path2mcscanx <- "~/path/to/MCScanX/"

library(genespace)
gpar <- init_genespace(
  wd = output_dir,
  genomeIDs = spec,
  ploidy = ploi,
  ignoreTheseGenomes = outgroup,
  minPepLen = 30,
  nCores = 48,
  maxOgPlaces = 24,
  orthofinderInBlk = FALSE
  )

gpar <- run_genespace(gsParam = gpar)

pangenome <- query_pangenes(
  gsParam,
  bed = NULL,
  refGenome = "AGB",
  transform = TRUE,
  showArrayMem = TRUE,
  showNSOrtho = TRUE,
  maxMem2Show = Inf,
  showUnPlacedPgs = TRUE
)

fwrite(pangenome,file="pangenome_AGB_synOG.txt",quote = F, sep = "\t", na = NA)
```

#### 3. Run the pipeline
The whole pipeline is organsed in several steps in jupty notebook. This requires that you have some knowledge of Python programming, however, this gives you a lot of freedom to modify the process as you like.


## Output
#### 1. process synOG
1) read in synOG table and adjust name
#### 2. perform drimm-synteny to call synteny blocks
1) process synOG table, call pairwise orthologs for given ref-target pair
2) DRIMM-Synteny was designed to handle complex genomes [23], for example genomes with duplications followed by diploidization. In the case of the comparison of S. taccada and Asteraceae diploid genomes, single, duplicated and triplicated synteny blocks in Asteraceae diploid genomes were reconstructed. For further analysis, the triplicated blocks (>= 5 anchor genes) that are presumably derived from the ancient Asteraceae genome triplication were used. ![image](https://github.com/user-attachments/assets/7762b1ab-490b-4ae0-889a-bb08a57829b7)

#### 3. process drimm raw blocks using processDrimm()
from drimm-synteny, we get raw blocks, and we need to clean and sort the raw blcoks to
1) split the all-in-one blocks (drimm-synteny output) into individual genomes
2) meanwhile divide the blocks into separate ones based on ratio (normally from genome duplication, but segmental duplications cannot be distinguished, and are also processed)
#### 4. with the results from step 3 (the final blocks), we can do:
1) ancestral genome reconstruction (multi-genome analysis)
   >a. with the outputs from step3, we can do ancestral chromosome reconstruction using IAGS
   >
   >b. different scenarios need to be considered based on the genome evolution history (duplication), and outgroup is very important
2) chracterize genome rearrangements, there are two approaches:
   >a. iags (fission, fusion), please follow iags guide; the output from step 3 can be used directly
   >
   >b. Grimm approach (fission, fusion, inversion, translocation)

## TODOs

## Citation
