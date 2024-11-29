# Asteraceae-synteny-phylogenomics
A pipeline for Asteraceae synteny phylogenomic analysis

### Installation

#### 1. set up an Linux version drimm-synteny <https://github.com/xjtu-omics/processDrimm/tree/master>
1) install mono <conda install mono -c conda-forge>; this is required to compile and run drimm-synteny
2) download the Program.cs, we used <mcs Program.cs -out:drimm-synteny> to compile a linux drimm-synteny; be noted that we can even modify the Program.cs to add a --help message
3) run drimm-synteny on linux: mono drimm-synteny [4 arguements]
