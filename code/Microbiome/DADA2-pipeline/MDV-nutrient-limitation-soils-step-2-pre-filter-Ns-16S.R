#Begin script 
#Load packages:
library(dada2); packageVersion("dada2") # the dada2 pipeline
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this
library(dplyr); packageVersion("dplyr") # for manipulating data

# Set path to project 
project.fp <- "/projects/adso6676/MDV_nutrient_limitation_soils/16s" # CHANGE ME to your project directory for any output you generate

# Set up names of sub directories
preprocess.fp <- file.path(project.fp, "01_preprocess")
demultiplex.fp <- file.path(preprocess.fp, "demultiplexed")
filtN.fp <- file.path(preprocess.fp, "filtN")
trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter")
table.fp <- file.path(project.fp, "03_tabletax")

# Set path for outputs from previous step in pipeline (i.e., demultiplex with idemp)
# Forward and reverse fastq filenames have format: 
fnFs <- sort(list.files(demultiplex.fp, pattern="R1_", full.names = TRUE))
fnRs <- sort(list.files(demultiplex.fp, pattern="R2_", full.names = TRUE))

################# Filter N ###########################
# Rename files - use gsub to get names in order!
R1_names <- gsub(paste0(demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_"), "", 
                 list.files(demultiplex.fp, pattern="R1_", full.names = TRUE))
print("See here R1 names")
print(R1_names)
file.rename(list.files(demultiplex.fp, pattern="R1_", full.names = TRUE), 
            paste0(demultiplex.fp, "/R1_", R1_names))

R2_names <- gsub(paste0(demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_"), "", 
                 list.files(demultiplex.fp, pattern="R2_", full.names = TRUE)) 
file.rename(list.files(demultiplex.fp, pattern="R2_", full.names = TRUE),
            paste0(demultiplex.fp, "/R2_", R2_names))

# Pre-filter to remove sequence reads with Ns
# Ambiguous bases will make it hard for cutadapt to find short primer sequences in the reads.
# To solve this problem, remove sequences with ambiguous bases (Ns)
# Name the N-filtered files to put them in filtN/ subdirectory
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs))
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))
print("HERE_print_fnFs_Rs")
print(fnFs.filtN)
print(fnRs.filtN)

# Filter Ns from reads and put them into the filtN directory
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, verbose = TRUE) 
# CHANGE multithread to FALSE on Windows (here and elsewhere in the program)