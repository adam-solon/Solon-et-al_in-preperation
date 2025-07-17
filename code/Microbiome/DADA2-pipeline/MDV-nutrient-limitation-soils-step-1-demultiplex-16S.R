# Begin script 
# Set up pathway to idemp (demultiplexing)
idemp <- "/projects/adso6676/software/anaconda/envs/ampliconEnv/idemp/idemp" #this is the path for the Rstudio server. CHANGE for other machines
system2(idemp) # Check that idemp is in your path and you can run shell commands from R

# Set path to lab data 
data.fp <- "/pl/active/schmidt-extrm/lab_sequences/2020_amplicon_miseq_summer/16S" #wherever your data is
# List all files in shared folder to check path
list.files(data.fp)

# Set path to project 
project.fp <- "/projects/adso6676/MDV_nutrient_limitation_soils/16s" # CHANGE ME to your project directory for any output you generate
barcode.fp <- "/projects/adso6676/MDV_nutrient_limitation_soils/barcodes"
list.files(barcode.fp)

# Set file paths for barcodes file, map file, and raw DNA sequences (fastqs)
barcode.fp <- file.path(barcode.fp, "may2023_16s_barcodes_modRQ.txt") # .txt file: barcode </t> sampleID
map.fp <- file.path(data.fp, "SampleSheetUsed (1).csv")
I1.fp <- file.path(data.fp, "Undetermined_S0_L001_I1_001.fastq.gz")
R1.fp <- file.path(data.fp, "Undetermined_S0_L001_R1_001.fastq.gz")
R2.fp <- file.path(data.fp, "Undetermined_S0_L001_R2_001.fastq.gz")

# Set up names of sub-directories
preprocess.fp <- file.path(project.fp, "01_preprocess")
demultiplex.fp <- file.path(preprocess.fp, "demultiplexed")
filtN.fp <- file.path(preprocess.fp, "filtN")
trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter")
table.fp <- file.path(project.fp, "03_tabletax")

################# Demultiplex #######################
# Call demultiplexing script - a pre-processing step before running DADA2
# Demultiplexing subsets out reads barcoded for samples from this specific project 
flags <- paste("-b", barcode.fp, "-I1", I1.fp, "-R1", R1.fp, "-R2", R2.fp, "-o", demultiplex.fp)
#run command
system2(idemp, args = flags)

# Look at output of demultiplexing
list.files(demultiplex.fp) # there are sample names.

# Clean up the output from idemp
# Change names of unassigned reads so they are not included in downstream processing
unassigned_1 <- paste0("mv", " ", demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_unsigned.fastq.gz",
                       " ", demultiplex.fp, "/Unassigned_reads1.fastq.gz")

unassigned_2 <- paste0("mv", " ", demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_unsigned.fastq.gz", 
                       " ", demultiplex.fp, "/Unassigned_reads2.fastq.gz")

#run command
system(unassigned_1)
system(unassigned_2)