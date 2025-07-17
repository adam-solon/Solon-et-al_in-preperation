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

# Set path for outputs from previous steps in pipeline (i.e., primer removal with cutadapt)
subF.fp <- file.path(filter.fp, "preprocessed_F") 
subR.fp <- file.path(filter.fp, "preprocessed_R") 
fastqFs <- sort(list.files(subF.fp, pattern="fastq.gz"))
fastqRs <- sort(list.files(subR.fp, pattern="fastq.gz"))
filtpathF <- file.path(subF.fp, "filtered") # files go into preprocessed_F/filtered/
filtpathR <- file.path(subR.fp, "filtered")

########################## Trim Reads ##############################
print("filterAndTrim")

#The maxEE and truncLen values, will most likely need to be change, especially after viewing your quality plots
#check the dada2 tutorial (not Fierer Lab) for details on these values and how to change them appropriately
#Options for maxEE that I have tried before: c(2,2), c(2,5), c(2,7) or Inf, wich means infinite.
filt_out <- filterAndTrim(fwd=file.path(subF.fp, fastqFs), filt=file.path(filtpathF, fastqFs),
                          rev=file.path(subR.fp, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                          truncLen=c(125,125), maxEE=Inf, truncQ=0, maxN=0, rm.phix=TRUE,
                          compress=TRUE, verbose=TRUE, multithread=FALSE)

# look at how many reads were kept
head(filt_out) 

# summary of samples in filt_out by percentage: between 0 and 95% remained, mean of 71 and median 79. with tl=125 and maxee=2
# changed tl=120 and now 0.35-95% remain, with med 88% and mean 78. I guess there are just a few crummy samples?
filt_out %>% 
  data.frame() %>% 
  mutate(Samples = rownames(.),
         percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything()) %>%
  summarise(min_remaining = paste0(round(min(percent_kept), 2), "%"), 
            median_remaining = paste0(round(median(percent_kept), 2), "%"),
            mean_remaining = paste0(round(mean(percent_kept), 2), "%"), 
            max_remaining = paste0(round(max(percent_kept), 2), "%"))
