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
# metadata=metadata[which(metadata$Mastic_product=="Mastic gum"),]

rownames(metadata)=metadata$Sample_name
counts_table2=counts_table[,which(colnames(counts_table)%in%metadata$Sample_name)]
dim(counts_table2)

y <- DGEList(counts_table2)
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

#keep genes that are present in at least 18 samples with at least 10 reads
keep.exprs <- filterByExpr(y, group=y$samples$Mastic_product, min.count=counts.cutoff, large.n=18, min.prop=1)
y <- y[keep.exprs,, keep.lib.sizes=TRUE]   #change to FALSE
dim(y)

#normalization between libraries
y <- calcNormFactors(y)
y$samples
plotMDS(y, col = as.numeric(y$samples$Timepoint))


#design matrix
condition <- factor(metadata$Mastic_product, levels = c("Placebo gum", "Mastic oil", "Mastic gum"))
levels(condition)=c("placebo_gum", "mastic_oil", "mastic_gum")
day = factor(metadata$Timepoint)
group=factor(paste(condition, day, sep = "_"))
table(metadata$Batch)
metadata$Batch=factor(metadata$Batch, levels=c("b","c","d"))
design_multi=model.matrix(~0+group+metadata$Batch)
colnames(design_multi)=c(levels(group),levels(metadata$Batch)[2:3])

y2 <- voom(y, design_multi, plot = T)
y2
eset=data.frame(y2)
# t2=which(metadata$Timepoint==2 & metadata$Mastic_product== "Mastic gum")
# eset[1,t2]

#Estimate the correlation between measurements made on the same subject
corfit= duplicateCorrelation(eset, design_multi, block=metadata$Patient_id)

#This inter-subject correlation is input into the linear model fit
corfit$consensus

#Create the linear model 
fit <- lmFit(y2,design_multi,block=metadata$Patient,correlation=corfit$consensus)

#make any comparisons between the different groups ,
cm <- makeContrasts(MasticGum.Day21vsMasticGum.Day0=(mastic_gum_3-mastic_gum_1),
                    MasticOil.Day21vsMasticOil.Day0=(mastic_oil_3-mastic_oil_1),
                    PlaceboGum.Day21vsPlaceboGum.Day0=(placebo_gum_3-placebo_gum_1),
                    MasticGum.Day7vsMasticGum.Day0=(mastic_gum_2-mastic_gum_1),
                    MasticOil.Day7vsMasticOil.Day0=(mastic_oil_2-mastic_oil_1),
                    PlaceboGum.Day7vsPlaceboGum.Day0=(placebo_gum_2-placebo_gum_1),
                    levels=design_multi)

#compute these contrasts and moderated t-tests
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
top.table <- topTable(fit2, coef= "MasticGum.Day21vsMasticGum.Day0", n=Inf)
write.table(top.table, file = "./limma_MasticGum.Day21vsMasticGum.Day0_results.tsv", sep ="\t", quote = FALSE, row.names = TRUE)
top.table <- topTable(fit2, coef= "MasticGum.Day7vsMasticGum.Day0", n=Inf)
write.table(top.table, file = "./limma_MasticGum.Day7vsMasticGum.Day0_results.tsv", sep ="\t", quote = FALSE, row.names = TRUE)

top.table <- topTable(fit2, coef= "MasticOil.Day21vsMasticOil.Day0", n=Inf)
write.table(top.table, file = "./limma_MasticOil.Day21vsMasticOil.Day0_results.tsv", sep ="\t", quote = FALSE, row.names = TRUE)
top.table <- topTable(fit2, coef= "MasticOil.Day7vsMasticOil.Day0", n=Inf)
write.table(top.table, file = "./limma_MasticOil.Day7vsMasticOil.Day0_results.tsv", sep ="\t", quote = FALSE, row.names = TRUE)

top.table <- topTable(fit2, coef= "PlaceboGum.Day21vsPlaceboGum.Day0", n=Inf)
write.table(top.table, file = "./limma_PlaceboGum.Day21vsPlaceboGum.Day0_results.tsv", sep ="\t", quote = FALSE, row.names = TRUE)
top.table <- topTable(fit2, coef= "PlaceboGum.Day7vsPlaceboGum.Day0", n=Inf)
write.table(top.table, file = "./limma_PlaceboGum.Day7vsPlaceboGum.Day0_results.tsv", sep ="\t", quote = FALSE, row.names = TRUE)

print(summary(decideTests(fit2)))

write.table(summary(decideTests(fit2)), file = "./limma_allComparisons_summary.tsv", sep ="\t", quote = FALSE)

