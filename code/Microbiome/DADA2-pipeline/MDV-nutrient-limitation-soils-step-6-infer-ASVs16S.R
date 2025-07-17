#Begin Script 
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

# Set path for outputs from previous steps in pipeline (i.e., trim reads with DADA2)
subF.fp <- file.path(filter.fp, "preprocessed_F") 
subR.fp <- file.path(filter.fp, "preprocessed_R") 
fastqFs <- sort(list.files(subF.fp, pattern="fastq.gz"))
fastqRs <- sort(list.files(subR.fp, pattern="fastq.gz"))
filtpathF <- file.path(subF.fp, "filtered") # files go into preprocessed_F/filtered/
filtpathR <- file.path(subR.fp, "filtered")

################################ Infer Amplicon Sequence Variants (ASV)###################################
# Housekeeping step - set up and verify the file names for the output:
# File parsing
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

# Sample names in order
sample.names <- substring(basename(filtFs), regexpr("_", basename(filtFs)) + 1) # doesn't drop fastq.gz
sample.names <- gsub(".fastq.gz", "", sample.names)
sample.namesR <- substring(basename(filtRs), regexpr("_", basename(filtRs)) + 1) # doesn't drop fastq.gz
sample.namesR <- gsub(".fastq.gz", "", sample.namesR)

# Double check
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#' #### Learn the error rates
set.seed(1603) # set seed to ensure that randomized steps are replicatable

# Learn forward error rates
errF <- learnErrors(filtFs, multithread=FALSE)

# Learn reverse error rates
errR <- learnErrors(filtRs, multithread=FALSE)


errF_plot <- plotErrors(errF, nominalQ=FALSE)
errR_plot <- plotErrors(errR, nominalQ=FALSE)

errF_plot #ok
errR_plot #ok

ggsave("16SerrF_plot.pdf", errF_plot, device="pdf")
ggsave("16serrR_plot.pdf", errR_plot, device="pdf")

# Dereplication, sequence inference, and merging of paired-end reads
# make lists to hold the loop output
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
ddF <- vector("list", length(sample.names))
names(ddF) <- sample.names
ddR <- vector("list", length(sample.names))
names(ddR) <- sample.names

# For each sample, get a list of merged and denoised sequences
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  # Dereplicate forward reads
  derepF <- derepFastq(filtFs[[sam]])
  # Infer sequences for forward reads
  dadaF <- dada(derepF, err=errF, multithread=TRUE)
  ddF[[sam]] <- dadaF
  # Dereplicate reverse reads
  derepR <- derepFastq(filtRs[[sam]])
  # Infer sequences for reverse reads
  dadaR <- dada(derepR, err=errR, multithread=TRUE)
  ddR[[sam]] <- dadaR
  # Merge reads together
  merger <- mergePairs(ddF[[sam]], derepF, ddR[[sam]], derepR)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# Save table as an r data object file
dir.create(table.fp)
saveRDS(seqtab, paste0(table.fp, "/seqtab.rds"))
saveRDS(mergers, paste0(table.fp, "/mergers.rds"))
