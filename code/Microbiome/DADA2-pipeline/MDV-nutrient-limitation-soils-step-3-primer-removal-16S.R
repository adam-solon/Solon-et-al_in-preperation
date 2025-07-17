#Install DADA2
#install.packages("BiocManager",repos = "http://cran.us.r-project.org")
#library("devtools")
#devtools::install_github("benjjneb/dada2", INSTALL_opts = '--no-lock')

#Load packages:
library(dada2); packageVersion("dada2") # the dada2 pipeline
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this
library(dplyr); packageVersion("dplyr") # for manipulating data

# Set up pathway to cutadapt (primer removal)
cutadapt <- "/projects/adso6676/software/anaconda/envs/ampliconEnv/bin/cutadapt" #this is the path for the Rstudio server. CHANGE for other machines
system2(cutadapt, args = "--version") # Check by running shell command from R

# Set path to project 
project.fp <- "/projects/adso6676/MDV_nutrient_limitation_soils/16s" # CHANGE ME to your project directory for any output you generate

# Set up names of sub directories 
preprocess.fp <- file.path(project.fp, "01_preprocess")
demultiplex.fp <- file.path(preprocess.fp, "demultiplexed")
filtN.fp <- file.path(preprocess.fp, "filtN")
trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter")
table.fp <- file.path(project.fp, "03_tabletax")

# Set path for outputs from previous step in pipeline (i.e., filter N with DADA2)
fnFs <- sort(list.files(demultiplex.fp, pattern="R1_", full.names = TRUE))
fnRs <- sort(list.files(demultiplex.fp, pattern="R2_", full.names = TRUE))
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs))
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))

#################### Primer Removal #######################################
#Set up the primer sequences to pass along to cutadapt
FWD <- "GTGYCAGCMGCCGCGGTAA"  ### this is 515f
REV <- "GGACTACNVGGGTWTCTAAT"  ### this is 806Br

# Write a function that creates a list of all orientations of the primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)# The Biostrings works w/ DNAString objects rather than character vectors
  #print(dna)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Save the primer orientations to pass to cutadapt
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
print(FWD.orients)

# Write a function that counts how many time primers appear in a sequence
print("Primer Hits")
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
print(fnFs.filtN[1])
print(fnFs.filtN[3])
print(fnFs.filtN[8])

#Before ruining cutadapt, we will look at primer detection for the first sample, as a check. 
#There may be some primers here, we will remove them below using CUTADAPT
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[3]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[3]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[3]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[3]]))
# 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[8]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[8]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[8]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[8]]))

# Remove primers with cutadapt and assess the output
# But first: Create directory to hold the output from cutadapt
if (!dir.exists(trimmed.fp)) dir.create(trimmed.fp)
print(fnFs)
fnFs.cut <- file.path(trimmed.fp, basename(fnFs))
fnRs.cut <- file.path(trimmed.fp, basename(fnRs))

# Save the reverse complements of the primers to variables
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Create the cutadapt flags ##
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 50")

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 50")
print(seq_along(fnFs)) #To see how many object are inside of the fnFs TO BE SURE that all of demultiplexed files are going to pass cutadaapt.

# Run Cutadapt
for (i in seq_along(fnFs)) {
  print("##############")
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#As a sanity check
print("Sanity check")
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[8]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[8]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[8]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[8]]))
