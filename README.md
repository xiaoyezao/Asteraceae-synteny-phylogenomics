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
