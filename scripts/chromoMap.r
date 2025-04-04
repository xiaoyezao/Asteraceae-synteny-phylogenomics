library(chromoMap)
## Set working directory
setwd("/Users/fengtao/Library/CloudStorage/OneDrive-WageningenUniversity&Research/project/Asteraceae_evolution/Asteraceae_phylogenomics/pipeline/")

## step 1. Load ALG color code
chr_colour_a <- read.table("sp03v1_chr_color_a.txt",header=T,sep="\t") # make sure the values in color column is surrounded by " "
chr_colour_a <- as.data.frame(chr_colour_a)

chr_colour_b <- read.table("sp03v1_chr_color_b.txt",header=T,sep="\t") # make sure the values in color column is surrounded by " "
chr_colour_b <- as.data.frame(chr_colour_b)

chr_colour_c <- read.table("sp03v1_chr_color_c.txt",header=T,sep="\t") # make sure the values in color column is surrounded by " "
chr_colour_c <- as.data.frame(chr_colour_c)

chr_colour <- read.table("sp03v1_chr_color.txt",header=T,sep="\t") # make sure the values in color column is surrounded by " "

#---------the data like this-------------
# drimm_chr group hex
# chr1  a1  #DF1159
# chr2  b2  #26AF67
# ...  ...
#----------------------

## step 2. Load chromosome information
setwd("/Users/fengtao/Library/CloudStorage/OneDrive-WageningenUniversity&Research/script/PhyloGenomics/AstGenEvo_paper/AGB_pipeline/workdir/sp46/6_chromPainting")
sp46_chr_info <- read.table("sp46.txt", header=T, sep="\t")
sp46_chr_info <- as.data.frame(sp46_chr_info)
#---------the data like this-------------
# chr haplotype start end ancestral_group
# chr_1 0  15  561 chr_15
# ...  ...
#----------------------


## step 3. Load segments information
setwd("/Users/fengtao/Library/CloudStorage/OneDrive-WageningenUniversity&Research/project/Asteraceae_evolution/2_analysis/10_ancestralGenome/genespace26/refSp03v1/sp14/6_chromPainting")
sp46_1_cooGene <- read.table("sp46_1_cooGene.txt", header=T, sep="\t")
sp46_1_cooGene <- as.data.frame(sp46_1_cooGene)
#---------the data like this-------------
# drimm_chr genome_chr start end
# chr_1 Alap.Chr01  0 4432
# ...  ...
#----------------------


## step 4. Plot
chromoMap(list(sp46_chr_info[, c("drimm_chr", "start", "end")]),
          list(sp46_1_cooGene[, c("chr", "chr","start", "end","ancestral_group")]),
          n_win.factor = 1,
          win.summary.display=T,
          segment_annotation=T,
          data_based_color_map=T,
          chr_color="grey",
          data_type="categorical",
          data_colors = list(chr_colour$hex),
          discrete.domain = list(chr_colour$drimm_chr),
          legend=T,
          lg_x=0,lg_y = 750,
          export.options=T,
          interactivity=F,
          scale.suffix = "",
          guides = T
          )
