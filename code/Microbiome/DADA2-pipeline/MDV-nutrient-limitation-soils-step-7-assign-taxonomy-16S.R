#DADA2
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

################### Remove Chimeras ########################
#call seqtab as st.all

st.all <- readRDS(paste0(table.fp, "/seqtab.rds"))
print("Here st.all")

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

print("porcentage of not chimeric")

saveRDS(seqtab.nochim, paste0(table.fp, "/seqtab.nochim.rds"))

print("Here after removing chemeras")


#print percentage of our sequences that were not chimeric.
100*sum(seqtab.nochim)/sum(st.all)

################### Assign Taxonomy ########################

tax <- assignTaxonomy(seqtab.nochim, 
                      "/pl/active/schmidt-extrm/databases/16S/silva_nr99_v138.1_train_set.fa.gz", 
                      tryRC = TRUE, multithread=TRUE) #set your own file path
#save results
saveRDS(seqtab.nochim, paste0(table.fp, "/seqtab_final.rds"))
saveRDS(tax, paste0(table.fp, "/tax_final.rds"))
write.table(seqtab.nochim,file = "seqtab_nochim.txt",sep="\t",quote=F,col.names=NA)

seqtab.nochim <- readRDS(paste0(table.fp, "/seqtab_final.rds"))
tax <- readRDS(paste0(table.fp, "/tax_final.rds"))

#setwd("/projects/adso6676/")