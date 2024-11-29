# Asteraceae-synteny-phylogenomics
A pipeline for Asteraceae synteny phylogenomic analysis

### Installation

#### 1. set up an Linux version drimm-synteny <https://github.com/xjtu-omics/processDrimm/tree/master>
1) install mono <conda install mono -c conda-forge>; this is required to compile and run drimm-synteny
2) download the Program.cs, we used <mcs Program.cs -out:drimm-synteny> to compile a linux drimm-synteny; be noted that we can even modify the Program.cs to add a --help message
3) run drimm-synteny on linux: mono drimm-synteny [4 arguements]

### Steps
#### 1. process synOG
1) read in synOG table and adjust name
#### 2. perform drimm-synteny to call synteny blocks
1) process synOG table, call pairwise orthologs for given ref-target pair
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
