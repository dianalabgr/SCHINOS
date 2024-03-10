#@@@@@@@@@@@@@ This script is used to normalize the data with CSS normalization  @@@@@@@@@@@@@@@@@@@@@@@@
# https://www.metagenomics.wiki/tools/16s/norm/css

library(metagenomeSeq); packageVersion("metagenomeSeq")

# other Packages
library("tidyverse"); packageVersion("tidyverse") 
library("coin"); packageVersion("coin") 
library('dplyr'); packageVersion("dplyr") 
library("pROC"); packageVersion("pROC") 
library("yaml"); packageVersion("yaml") 

setwd("/mnt/raid1/armen/Schinos/META-ANALYSIS/ALL/Fungi/")
getwd()

#import meta and data
meta = read.table(file             = 'RNA_metadata_ALL.tab', 
                  header           = TRUE, 
                  sep              = "\t", 
                  row.names=1,
                  stringsAsFactors = FALSE)
#meta[1:5, ]

#here, the abundance data of genera is used
data = read.table(file             = "genus.tab",
                  header           = TRUE,
                  sep              = "\t", 
                  row.names=1,
                  stringsAsFactors = FALSE)

#colnames(data)=rownames(meta)

data<- t(data)

dim(data)
head(data)

# convert into package format
metaSeqObject = newMRexperiment(data)
# CSS normalization
# cumNormStat instead of cumNormStatFast
metaSeqObject_CSS  = cumNorm(metaSeqObject, p=cumNormStat(metaSeqObject))
# convert CSS normalized data into data.frame-formatted table 
data_CSS = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
data_CSS<- t(data_CSS)
# Save the dataframe as a tab-delimited file
write.table(data_CSS, file = "genus_CSS.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
