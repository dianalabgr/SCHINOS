
# BiocManager::install(version = "3.16")
library("BiocManager")
# BiocManager::install("edgeR")
library(edgeR)
#BiocManager::install("tximportData")


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
keep_1=which(metadata$BoP>30 & metadata$Timepoint==1)
keep_3=which(metadata$BoP<=30 & metadata$Timepoint==3)

metadata_keep=metadata[c(keep_1,keep_3),]
dim(metadata_keep)
unique(metadata_keep$Patient_id)
metadata_keep=metadata_keep[-which(metadata_keep$Patient_id%in%c("S5","S11")),]

rownames(metadata_keep)=metadata_keep$Sample_name
counts_table2=counts_table[,which(colnames(counts_table)%in%metadata_keep$Sample_name)]
dim(counts_table2)

##################################################################
########### EdgeR with paired comparison ##########################

y <- DGEList(counts_table2)
y$samples=merge(y$samples,metadata_keep, by="row.names")
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
table(y$samples$Timepoint)
#keep genes that are present in at least 18 samples with at least 10 reads
keep.exprs <- filterByExpr(y, group=y$samples$Timepoint, min.count=counts.cutoff, large.n=18, min.prop=1)
y <- y[keep.exprs,, keep.lib.sizes=TRUE]   #change to FALSE
dim(y)


#normalization between libraries
y <- calcNormFactors(y)
metadata = y$samples
table(metadata$Timepoint, metadata$Mastic_product)

#design matrix
#page 43 in edgeR manual
subjectid=factor(metadata$Patient_id)
day = factor(metadata$Timepoint)
batch=metadata$Batch
batch=factor(batch, levels=c("b","c","d"))
lib.size=metadata$lib.size

dim(metadata)
design <- model.matrix(~subjectid+day)

#Dispersion estimation
y <- estimateDisp(y, design, robust=TRUE)
# plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

#RNAs different in day 21 from 0 in placebo
qlf <- glmQLFTest(fit, coef="day3")
topTags(qlf)
is.de <- decideTestsDGE(qlf)
summary(is.de)
write.table(summary(is.de), file = "./edgeR_BoP_t3.vs.t1_summary.tsv", sep ="\t", quote = FALSE)
result_table = topTags(qlf, n = Inf)$table
write.table(result_table, file = "./edgeR_BoP_t3.vs.t1_DE.tsv", sep ="\t", quote = FALSE, row.names = TRUE)


####
#Batch effect 

design <- model.matrix(~day+batch)

#Dispersion estimation
y <- estimateDisp(y, design, robust=TRUE)
# plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

#RNAs different in day 21 from 0 in placebo
qlf <- glmQLFTest(fit, coef="day3")
topTags(qlf)
is.de <- decideTestsDGE(qlf)
summary(is.de)
write.table(summary(is.de), file = "./edgeR_S2_BoP_t3.vs.t1_summary.tsv", sep ="\t", quote = FALSE)
result_table = topTags(qlf, n = Inf)$table
write.table(result_table, file = "./edgeR_S2_BoP_t3.vs.t1_DE.tsv", sep ="\t", quote = FALSE, row.names = TRUE)

#################################################
###################### Limma ####################

y <- DGEList(counts_table2)
y$samples=merge(y$samples,metadata_keep, by="row.names")
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
keep.exprs <- filterByExpr(y, group=y$samples$Timepoint, min.count=counts.cutoff, large.n=18, min.prop=1)
y <- y[keep.exprs, keep.lib.sizes=TRUE]   #change to FALSE
dim(y)

#normalization between libraries
y <- calcNormFactors(y)
y$samples
plotMDS(y, col = as.numeric(y$samples$Timepoint))


#design matrix
day = factor(metadata_keep$Timepoint)
levels(day)=c("day0", "day21")
table(metadata_keep$Batch)
metadata_keep$Batch=factor(metadata_keep$Batch, levels=c("b","c","d"))

design_multi=model.matrix(~0+day+metadata_keep$Batch)
colnames(design_multi)=c(levels(day),levels(metadata_keep$Batch)[2:3])

y2 <- voom(y, design_multi, plot = T)
y2
eset=data.frame(y2)
# t2=which(metadata$Timepoint==2 & metadata$Mastic_product== "Mastic gum")
# eset[1,t2]

#Estimate the correlation between measurements made on the same subject
corfit= duplicateCorrelation(eset, design_multi, block=metadata_keep$Patient_id)

#This inter-subject correlation is input into the linear model fit
corfit$consensus

#Create the linear model 
fit <- lmFit(y2,design_multi,block=metadata$Patient,correlation=corfit$consensus)

#make any comparisons between the different groups ,
cm <- makeContrasts(day21_vs_day0=(day21-day0),
                    levels=design_multi)

#compute these contrasts and moderated t-tests
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
top.table <- topTable(fit2, coef= "day21_vs_day0", n=Inf)
print(summary(decideTests(fit2)))
write.table(summary(decideTests(fit2)), file = "./limma_BoP_t3.vs.t1_summary.tsv", sep ="\t", quote = FALSE)

write.table(top.table, file = "./limma_BoP_t3vst1_results.tsv", sep ="\t", quote = FALSE, row.names = TRUE)

        