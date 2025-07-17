#Install DADA2
#install.packages("BiocManager",repos = "http://cran.us.r-project.org")
#library("devtools")
#devtools::install_github("benjjneb/dada2", INSTALL_opts = '--no-lock')

#Load packages:
library(dada2); packageVersion("dada2") # the dada2 pipeline
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this
library(dplyr); packageVersion("dplyr") # for manipulating data
library(tidyr); packageVersion("tidyr") # for graphs
library(Hmisc); packageVersion("Hmisc") # for graphs
library(ggplot2); packageVersion("ggplot2") # for graphs
library(plotly); packageVersion("plotly") # for interactive graphs 

# Set path to project 
project.fp <- "/projects/adso6676/MDV_nutrient_limitation_soils/16s" # CHANGE ME to your project directory for any output you generate

# Set up names of sub directories 
preprocess.fp <- file.path(project.fp, "01_preprocess")
demultiplex.fp <- file.path(preprocess.fp, "demultiplexed")
filtN.fp <- file.path(preprocess.fp, "filtN")
trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter")
table.fp <- file.path(project.fp, "03_tabletax")

# Set path for outputs from previous step in pipeline (i.e., primer removal with cutadapt)
fnFs <- sort(list.files(demultiplex.fp, pattern="R1_", full.names = TRUE))
fnRs <- sort(list.files(demultiplex.fp, pattern="R2_", full.names = TRUE))
fnFs.cut <- file.path(trimmed.fp, basename(fnFs))
fnRs.cut <- file.path(trimmed.fp, basename(fnRs))

###################### Assess Quality Scores ######################################
#Put filtered reads into separate sub-directories for big data workflow
#Here new files will appear in the new 02_filter folder
dir.create(filter.fp)
subF.fp <- file.path(filter.fp, "preprocessed_F") 
subR.fp <- file.path(filter.fp, "preprocessed_R") 
dir.create(subF.fp)
dir.create(subR.fp)

# Move R1 and R2 from trimmed to separate forward/reverse sub-directories
fnFs.Q <- file.path(subF.fp,  basename(fnFs)) 
fnRs.Q <- file.path(subR.fp,  basename(fnRs))
file.rename(from = fnFs.cut, to = fnFs.Q)
file.rename(from = fnRs.cut, to = fnRs.Q)

# File parsing; create file names and make sure that forward and reverse files match
filtpathF <- file.path(subF.fp, "filtered") # files go into preprocessed_F/filtered/
filtpathR <- file.path(subR.fp, "filtered") # files go into preprocessed_F/filtered/
fastqFs <- sort(list.files(subF.fp, pattern="fastq.gz"))
fastqRs <- sort(list.files(subR.fp, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#Run quality plots
#Quality plots for a selection of samples in order to help you determine what # of bp to cut off at (160 was mine) 
#If the number of samples is 20 or less, plot them all, otherwise, just plot 20 randomly selected samples
if( length(fastqFs) <= 20) {
  fwd_qual_plots <- plotQualityProfile(paste0(subF.fp, "/", fastqFs))
  rev_qual_plots <- plotQualityProfile(paste0(subR.fp, "/", fastqRs))
} else {
  rand_samples <- sample(size = 20, 1:length(fastqFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(paste0(subF.fp, "/", fastqFs[rand_samples]))
  rev_qual_plots <- plotQualityProfile(paste0(subR.fp, "/", fastqRs[rand_samples]))
}
par(mar=rep(4,4))

fwd_qual_plots #look ok out to probably 120 or farther
rev_qual_plots 

#write plots to my project path bc it won't be visualized.
saveRDS(fwd_qual_plots, paste0(filter.fp, "/fwd_qual_plots16S.rds"))
saveRDS(rev_qual_plots, paste0(filter.fp, "/rev_qual_plots16S.rds"))

ggsave("Svalbard_2023_18S/fwd_qual_plots16S.pdf", fwd_qual_plots, device="pdf")
ggsave("Svalbard_2023_18S/rev_qual_plots16S.pdf", rev_qual_plots, device="pdf")
