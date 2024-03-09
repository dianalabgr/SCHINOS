library(dplyr)
# BiocManager::install(version = "3.16")
library("BiocManager")
# BiocManager::install("edgeR")
library(edgeR)
#BiocManager::install("tximportData")
library(tidyverse)


dir <- "/mnt/raid1/philip/schinos_seq/data/mrna/mapped_to_genes"
files = list.files(dir)[grepl("ReadsPerGene", list.files(dir))]
names(files) <- gsub("_S.*ReadsPerGene.*","",files)
files_paths=file.path(dir,files)
names(files_paths) <- names(files)

counts_table = sapply(files_paths, function(i)
{
  #the 10th column has the counts
  x = read.csv(i, sep="\t", stringsAsFactors = FALSE, header = FALSE, skip=4)[, 4]
})
rownames(counts_table) = rownames(read.csv(files_paths[1], sep="\t", stringsAsFactors = FALSE, header = FALSE, row.names = 1, skip=4))

metadata_file="/mnt/raid1b/philip/schinos_seq/metadata/SCHINOS_Sequenced_samples_metadata.tsv"
metadata = read.csv(metadata_file, sep="\t", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$Sample_name
metadata$Timepoint=factor(metadata$Timepoint)
metadata$Mastic_product=relevel(factor(metadata$Mastic_product),"Placebo gum")
metadata$Batch=factor(metadata$Batch)
metadata = metadata %>% filter(Discarded==0,na.rm=TRUE)
metadata = metadata %>% filter(str_detect(Sample_name,"MRN"),na.rm=TRUE)

#change order of metadata rows to match files list
metadata=metadata[match(names(files_paths),metadata$Sample_name),]
metadata=metadata[which(metadata$Sample_type_abbreviation!='TI'),]
#metadata=metadata[which(metadata$Timepoint!='2'),]

rownames(metadata)=metadata$Sample_name
counts_table2=counts_table[,which(colnames(counts_table)%in%metadata$Sample_name)]
dim(counts_table2)

y <- DGEList(counts_table2, group = metadata$Mastic_product)
y$samples=merge(y$samples,metadata, by="row.names")
y$samples
counts_table[which(rownames(counts_table)=="ENSG00000067167"),]
dim(y)

#choose minimum cutoff
counts.cutoff = 10

## Filtering and normalisation
libsizes = y$samples$lib.size
#keep samples with depth >=1 million
y = y[,which(y$samples$lib.size>1000000)]   
dim(y)
table(y$samples$group)
#keep genes that are present in at least 18 samples with at least 10 reads
keep.exprs <- filterByExpr(y, group=y$samples$group, min.count=counts.cutoff, large.n=18, min.prop=1)
y <- y[keep.exprs,, keep.lib.sizes=TRUE]   #change to FALSE
dim(y)


#normalization between libraries
y <- calcNormFactors(y)
metadata = y$samples
table(metadata$Timepoint, metadata$Mastic_product)

#design matrix
#page 43 in edgeR manual

######################################################################
########## Comparisons day21 vs day 0 
dim(metadata)
metadata_day21=metadata[which(metadata$Timepoint!='2'),]
subjectid=factor(metadata_day21$Patient_id)
day = factor(metadata_day21$Timepoint)
batch=metadata_day21$Batch
lib.size=metadata_day21$lib.size
condition <- factor(metadata_day21$Mastic_product, levels = c("Placebo gum", "Mastic oil", "Mastic gum"))
levels(condition)=c("placebo_gum", "mastic_oil", "mastic_gum")
group=factor(paste(condition, day, sep = "_"))
# List the current sample names
sample_names = colnames(y)
y_day21=y[,sample_names %in% metadata_day21$Row.names]

design <- model.matrix(~subjectid)
Placebo.Day21 = condition=="placebo_gum" & day=="3"
MasticOil.Day21 = condition=="mastic_oil" & day=="3"
MasticGum.Day21 = condition=="mastic_gum" & day=="3"

#design=cbind(design,Placebo.Day21,MasticOil.Day21,MasticGum.Day21)
design=cbind(design,Placebo.Day21,MasticOil.Day21,MasticGum.Day21)
#Dispersion estimation
y_day21 <- estimateDisp(y_day21, design, robust=TRUE)
# plotBCV(y)
fit <- glmQLFit(y_day21, design, robust=TRUE)
plotQLDisp(fit)

#RNAs different in day 21 from 0 in placebo
qlf <- glmQLFTest(fit, coef="Placebo.Day21")
topTags(qlf)
is.de <- decideTestsDGE(qlf)
summary(is.de)
write.table(summary(is.de), file = "./edgeR_mg_t3.vs.t1_summary_placebo.tsv", sep ="\t", quote = FALSE)
result_table = topTags(qlf, n = Inf)$table
write.table(result_table, file = "./edgeR_mg_t3.vs.t1_DE_placebo.tsv", sep ="\t", quote = FALSE, row.names = TRUE)


#RNAs different in day 21 from 0 in mastic gum
qlf <- glmQLFTest(fit, coef="MasticGum.Day21")
topTags(qlf)
is.de <- decideTestsDGE(qlf)
summary(is.de)
write.table(summary(is.de), file = "./edgeR_mg_t3.vs.t1_summary_masticGum.tsv", sep ="\t", quote = FALSE)
result_table = topTags(qlf, n = Inf)$table
write.table(result_table, file = "./edgeR_mg_t3.vs.t1_DE_masticGum.tsv", sep ="\t", quote = FALSE, row.names = TRUE)



#RNAs different in day 21 from 0 in mastic oil
qlf <- glmQLFTest(fit, coef="MasticOil.Day21")
topTags(qlf)
is.de <- decideTestsDGE(qlf)
summary(is.de)
write.table(summary(is.de), file = "./edgeR_mg_t3.vs.t1_summary_masticOil.tsv", sep ="\t", quote = FALSE)
result_table = topTags(qlf, n = Inf)$table
write.table(result_table, file = "./edgeR_mg_t3.vs.t1_DE_masticOil.tsv", sep ="\t", quote = FALSE, row.names = TRUE)


######################################################################
########## Comparisons day21 vs day 0 
dim(metadata)
metadata_day7=metadata[which(metadata$Timepoint!='3'),]
subjectid=factor(metadata_day7$Patient_id)
day = factor(metadata_day7$Timepoint)
batch=metadata_day7$Batch
lib.size=metadata_day7$lib.size
condition <- factor(metadata_day7$Mastic_product, levels = c("Placebo gum", "Mastic oil", "Mastic gum"))
levels(condition)=c("placebo_gum", "mastic_oil", "mastic_gum")
group=factor(paste(condition, day, sep = "_"))
# List the current sample names
sample_names = colnames(y)
y_day7=y[,sample_names %in% metadata_day7$Row.names]

design <- model.matrix(~subjectid)
Placebo.Day7 = condition=="placebo_gum" & day=="2"
MasticOil.Day7 = condition=="mastic_oil" & day=="2"
MasticGum.Day7 = condition=="mastic_gum" & day=="2"

design=cbind(design,Placebo.Day7,MasticOil.Day7,MasticGum.Day7)
#Dispersion estimation
y_day7 <- estimateDisp(y_day7, design, robust=TRUE)
# plotBCV(y)
fit <- glmQLFit(y_day7, design, robust=TRUE)
plotQLDisp(fit)



#RNAs different in day 7 from 0 in placebo
qlf <- glmQLFTest(fit, coef="Placebo.Day7")
topTags(qlf)
is.de <- decideTestsDGE(qlf)
summary(is.de)
write.table(summary(is.de), file = "./edgeR_mg_t2.vs.t1_summary_placebo.tsv", sep ="\t", quote = FALSE)
result_table = topTags(qlf, n = Inf)$table
write.table(result_table, file = "./edgeR_mg_t2.vs.t1_DE_placebo.tsv", sep ="\t", quote = FALSE, row.names = TRUE)


#RNAs different in day 7 from 0 in mastic gum
qlf <- glmQLFTest(fit, coef="MasticGum.Day7")
topTags(qlf)
is.de <- decideTestsDGE(qlf)
summary(is.de)
write.table(summary(is.de), file = "./edgeR_mg_t2.vs.t1_summary_masticGum.tsv", sep ="\t", quote = FALSE)
result_table = topTags(qlf, n = Inf)$table
write.table(result_table, file = "./edgeR_mg_t2.vs.t1_DE_masticGum.tsv", sep ="\t", quote = FALSE, row.names = TRUE)



#RNAs different in day 7 from 0 in mastic oil
qlf <- glmQLFTest(fit, coef="MasticOil.Day7")
topTags(qlf)
is.de <- decideTestsDGE(qlf)
summary(is.de)
write.table(summary(is.de), file = "./edgeR_mg_t2.vs.t1_summary_masticOil.tsv", sep ="\t", quote = FALSE)
result_table = topTags(qlf, n = Inf)$table
write.table(result_table, file = "./edgeR_mg_t2.vs.t1_DE_masticOil.tsv", sep ="\t", quote = FALSE, row.names = TRUE)



######################################################################
########## Comparisons day21 vs day 0 with batch effect
dim(metadata)
metadata_day21=metadata[which(metadata$Timepoint!='2'),]
subjectid=factor(metadata_day21$Patient_id)
day = factor(metadata_day21$Timepoint)
batch=metadata_day21$Batch
batch=factor(batch, levels=c("b","c", "d"))
lib.size=metadata_day21$lib.size
condition <- factor(metadata_day21$Mastic_product, levels = c("Placebo gum", "Mastic oil", "Mastic gum"))
levels(condition)=c("placebo_gum", "mastic_oil", "mastic_gum")
group=factor(paste(condition, day, sep = "_"))
# List the current sample names
sample_names = colnames(y)
y_day21=y[,sample_names %in% metadata_day21$Row.names]

design <- model.matrix(~subjectid+batch)
Placebo.Day21 = condition=="placebo_gum" & day=="3"
MasticOil.Day21 = condition=="mastic_oil" & day=="3"
MasticGum.Day21 = condition=="mastic_gum" & day=="3"

#design=cbind(design,Placebo.Day21,MasticOil.Day21,MasticGum.Day21)
design=cbind(design,Placebo.Day21,MasticOil.Day21,MasticGum.Day21)
#Dispersion estimation
y_day21 <- estimateDisp(y_day21, design, robust=TRUE)
# plotBCV(y)
fit <- glmQLFit(y_day21, design, robust=TRUE)
plotQLDisp(fit)


######################################################################
########## Comparisons day21 vs day 0 with batch effect 
dim(metadata)
metadata_day7=metadata[which(metadata$Timepoint!='3'),]
subjectid=factor(metadata_day7$Patient_id)
day = factor(metadata_day7$Timepoint)
batch=metadata_day7$Batch
batch = factor(batch, levels=c("b", "c", "d"))
lib.size=metadata_day7$lib.size
condition <- factor(metadata_day7$Mastic_product, levels = c("Placebo gum", "Mastic oil", "Mastic gum"))
levels(condition)=c("placebo_gum", "mastic_oil", "mastic_gum")
group=factor(paste(condition, day, sep = "_"))
# List the current sample names
sample_names = colnames(y)
y_day7=y[,sample_names %in% metadata_day7$Row.names]

design <- model.matrix(~subjectid+batch)
Placebo.Day7 = condition=="placebo_gum" & day=="2"
MasticOil.Day7 = condition=="mastic_oil" & day=="2"
MasticGum.Day7 = condition=="mastic_gum" & day=="2"

design=cbind(design,Placebo.Day7,MasticOil.Day7,MasticGum.Day7)
#Dispersion estimation
y_day7 <- estimateDisp(y_day7, design, robust=TRUE)
# plotBCV(y)
fit <- glmQLFit(y_day7, design, robust=TRUE)
plotQLDisp(fit)

        