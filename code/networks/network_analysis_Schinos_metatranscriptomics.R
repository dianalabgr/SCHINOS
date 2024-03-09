library(NetCoMi)
library(MetBrewer)
library(data.table)
library(wesanderson)
library(dplyr)
`%nin%`=Negate(`%in%`)

#Interesting bacteria from bibliografy
interesting_bacteria=c("Aggregatibacter actinomycetemcomitans", "Campylobacter rectus", "Eikenella corrodens",
                       "Filifactor alocis", "Fusobacterium nucleatum", "Parvimonas micra", 
                       "Porphyromonas gingivalis", "Prevotella intermedia", "Streptococcus intermedius",
                       "Tannerella forsythia", "Treponema denticola")


#Upload metadata and keep only metagenomcs metadata 
metadata=read.table(file="/mnt/raid1b/philip/schinos_seq/metadata/SCHINOS_Sequenced_samples_metadata.tsv", sep="\t", header=T)
metadata=metadata[-which(metadata$Discarded==1),]

#Keep the subjects with BOP > 30 in timepoint 1 and BoP <= 3- in Timepoint3 
metadata_timepoint1=metadata[which(metadata$Timepoint==1 & metadata$BoP > 30),]
metadata_timepoint3=metadata[which(metadata$Timepoint==3 & metadata$BoP <= 30),]
keep_subjects = intersect(unique(metadata_timepoint1$Patient_id),unique(metadata_timepoint3$Patient_id))
#Remove also G12 and G7 for multiple reasons
keep_subjects2= keep_subjects[which(keep_subjects%nin%c("G7","G12"))]
#Keep the metadata only for subjects in keep_subjects2 
metadata_timepoint1=metadata_timepoint1[which(metadata_timepoint1$Patient_id%in%keep_subjects2),]
metadata_timepoint3=metadata_timepoint3[which(metadata_timepoint3$Patient_id%in%keep_subjects2),]


abundances=read.delim(file="/mnt/raid1b/philip/schinos_seq/results/metatranscriptomics/abundances/metatranscriptomics_abundances.tsv")
column_sums=c()
column_sums_decimal=c()
# Compute the sum of each row for each column
colnames(abundances)[18]
for (i in seq(2,142,1)) {
  int=round(sum(abundances[,i]))
  column_sums=c(column_sums,int)
  column_sums_decimal=c(column_sums_decimal,sum(abundances[,i]))
}
column_sums
#column_sums <- apply(abundances_total, 2, sum)
#col_sums <- colSums(abundances_total)
write.csv(data.frame("samples"=colnames(abundances[2:142]), "agamemnon_pairs"=column_sums, 
                     "agamemnon_pairs_decimal"=column_sums_decimal), file="total_agamemnon_PairsCounts_metatranscriptomics.csv")


#Lets start with the SB tissue in timepoint 1 
abundances_SB_timepoint1=abundances[,which(colnames(abundances)%in%metadata_timepoint1[which(metadata_timepoint1$Sample_type_abbreviation=="SB"),"Sample_name"])]
rownames(abundances_SB_timepoint1)=abundances$Name
metadata_SB = metadata_timepoint1[which(metadata_timepoint1$Sample_type_abbreviation=="SB"),]

#Find the relative frequency of microbiomes in each group
abundances_SB_timepoint1_2= apply(abundances_SB_timepoint1,2,function(x) x/sum(x))
rownames(abundances_SB_timepoint1_2)=rownames(abundances_SB_timepoint1)
abundances_SB_timepoint1_3=abundances_SB_timepoint1_2
abundances_SB_timepoint1_3= cbind(abundances_SB_timepoint1_2,median = unlist(apply(abundances_SB_timepoint1_2, 1, median)))
abundances_SB_timepoint1_3 =  data.frame(cbind(abundances_SB_timepoint1_3,mean = unlist(apply(abundances_SB_timepoint1_2, 1, mean))))

abundances_SB_timepoint1_3[which(rownames(abundances_SB_timepoint1_3)=="Porphyromonas gingivalis"),]
abundances_SB_timepoint1_3[which(rownames(abundances_SB_timepoint1_3)=="Tannerella forsythia"),]

#Find which microbiomes we will keep
keep_SB_timepoint1_median = rownames(abundances_SB_timepoint1_3[order(abundances_SB_timepoint1_3$median, decreasing = TRUE), ])[1:100]
keep_SB_timepoint1_mean = rownames(abundances_SB_timepoint1_3[order(abundances_SB_timepoint1_3$mean, decreasing = TRUE), ])[1:100]

#Keep the 100 microbiomes with the highest median relative frequency in each group 
abundances_SB_timepoint1_4=abundances_SB_timepoint1[which(rownames(abundances_SB_timepoint1)%in%keep_SB_timepoint1_median), ]

#Create a network using Spieceasi with positive and negative connections
#samples as rows and microbiomes as columns
net_spieceasi_SB_timepoint1_allConnections <- netConstruct(t(abundances_SB_timepoint1_4),
                                                           dataType = "counts",
                                                           filtTax = c("none"),
                                                           measure = "spieceasi", 
                                                           measurePar=list(method = "mb", nlambda = 10, 
                                                                           lambda.min.ratio = 0.1),
                                                           jointPrepro = NULL, 
                                                           filtSamp = "none", filtSampPar = NULL,  zeroMethod = "none",
                                                           zeroPar = NULL, normMethod = "none", normPar = NULL, sparsMethod = "none",
                                                           cores = 20,  dissFunc = "signed",  weighted = TRUE, verbose = 2, 
                                                           seed = 123456)


#Analyze the network
props_pears_SB_timepoint1_allConnections <- netAnalyze(net_spieceasi_SB_timepoint1_allConnections, 
                                                       centrLCC = FALSE,
                                                       avDissIgnoreInf = TRUE,
                                                       sPathNorm = FALSE,
                                                       clustMethod = "cluster_fast_greedy",
                                                       hubPar = c("between", "degree"),
                                                       hubQuant = 0.85,
                                                       lnormFit = TRUE,
                                                       normDeg = TRUE,
                                                       normBetw = TRUE)

summary(props_pears_SB_timepoint1_allConnections)
summary(props_pears_SB_timepoint1_allConnections)
write.table(summary(props_pears_SB_timepoint1_allConnections)[1],  file("microNetProps_summary_SB_timepoint1.txt"))
write.table(summary(props_pears_SB_timepoint1_allConnections)[4],  file("microNetProps_hubs_SB_timepoint1.txt"))

feature_nodes_periodo = props_pears_SB_timepoint1_allConnections$clustering$clust1
write.csv(feature_nodes_periodo,  file("microNetProps_clusters_SB_timepoint1.csv"))
interesting_feautre_nodes = feature_nodes_periodo[which(names(feature_nodes_periodo)%in%interesting_bacteria)]
write.csv(interesting_feautre_nodes,  file("microNetProps_interesting_clusters_SB_timepoint1.csv"))

table(feature_nodes_periodo)

#Define colors used
colors_manual_peri=c(wes_palette("GrandBudapest2")[2],"#FFFF00", wes_palette("Darjeeling2")[4], 
                     wes_palette("Darjeeling2")[2],  "#FF0000",wes_palette("Darjeeling1")[4],  
                     wes_palette("Darjeeling1")[3], "#00FF00",wes_palette("FantasticFox1")[2], 
                     "#8A2BE2", wes_palette("Moonrise2")[3], wes_palette("Moonrise3")[1], wes_palette("IsleofDogs1")[5], 
                     wes_palette("Moonrise3")[2], wes_palette("Darjeeling1")[1], wes_palette("IsleofDogs1")[1], 
                     wes_palette("Moonrise3"),wes_palette("Darjeeling1"))
#Plot the network
png("SB_timepoint1_AllConnections_MedianAbundance100Manual_labda01_nlabda30_NoPulsar.png", width = 10, height = 10, units = "in", res=600)
plot(props_pears_SB_timepoint1_allConnections, 
     nodeColor = "feature", 
     featVecCol = feature_nodes_periodo,
     colorVec = colors_manual_peri,
     nodeSize = "CSS",
     cexNodes = 2,
     nodeSizeSpread = 2, 
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     repulsion = 1.1, 
     rmSingles = "none",
     nodeFilter = "none", 
     nodeFilterPar = 1, 
     nodeTransp = 10, 
     hubTransp = 10, 
     highlightHubs = T,
     cexLabels=4,
     showTitle = TRUE,
     title1 = "",
     cexTitle = 10)
dev.off()


#Lets continue with the SB tissue in timepoint 3
abundances_SB_timepoint3=abundances[,which(colnames(abundances)%in%metadata_timepoint3[which(metadata_timepoint3$Sample_type_abbreviation=="SB"),"Sample_name"])]
rownames(abundances_SB_timepoint3)=abundances$Name
metadata_SB = metadata_timepoint3[which(metadata_timepoint3$Sample_type_abbreviation=="SB"),]

#Find the relative frequency of microbiomes in each group
abundances_SB_timepoint3_2= apply(abundances_SB_timepoint3,2,function(x) x/sum(x))
rownames(abundances_SB_timepoint3_2)=rownames(abundances_SB_timepoint3)
abundances_SB_timepoint3_3=abundances_SB_timepoint3_2
abundances_SB_timepoint3_3= cbind(abundances_SB_timepoint3_2,median = unlist(apply(abundances_SB_timepoint3_2, 1, median)))
abundances_SB_timepoint3_3 =  data.frame(cbind(abundances_SB_timepoint3_3,mean = unlist(apply(abundances_SB_timepoint3_2, 1, mean))))

abundances_SB_timepoint3_3[which(rownames(abundances_SB_timepoint3_3)=="Porphyromonas gingivalis"),]
abundances_SB_timepoint3_3[which(rownames(abundances_SB_timepoint3_3)=="Edwardsiella ictaluri"),]

#Find which microbiomes we will keep
keep_SB_timepoint3_median = rownames(abundances_SB_timepoint3_3[order(abundances_SB_timepoint3_3$median, decreasing = TRUE), ])[1:100]
keep_SB_timepoint3_mean = rownames(abundances_SB_timepoint3_3[order(abundances_SB_timepoint3_3$mean, decreasing = TRUE), ])[1:100]

#Keep the 100 microbiomes with the highest median relative frequency in each group 
abundances_SB_timepoint3_4=abundances_SB_timepoint3[which(rownames(abundances_SB_timepoint3)%in%keep_SB_timepoint3_median), ]

#Create a network using Spieceasi with positive and negative connections
#samples as rows and microbiomes as columns
net_spieceasi_SB_timepoint3_allConnections <- netConstruct(t(abundances_SB_timepoint3_4),
                                                           dataType = "counts",
                                                           filtTax = c("none"),
                                                           measure = "spieceasi", 
                                                           measurePar=list(method = "mb", nlambda = 10, 
                                                                           lambda.min.ratio = 0.1),
                                                           jointPrepro = NULL, 
                                                           filtSamp = "none", filtSampPar = NULL,  zeroMethod = "none",
                                                           zeroPar = NULL, normMethod = "none", normPar = NULL, sparsMethod = "none",
                                                           cores = 20,  dissFunc = "signed",  weighted = TRUE, verbose = 2, 
                                                           seed = 123456)


#Analyze the network
props_pears_SB_timepoint3_allConnections <- netAnalyze(net_spieceasi_SB_timepoint3_allConnections, 
                                                       centrLCC = FALSE,
                                                       avDissIgnoreInf = TRUE,
                                                       sPathNorm = FALSE,
                                                       clustMethod = "cluster_fast_greedy",
                                                       hubPar = c("between", "degree"),
                                                       hubQuant = 0.85,
                                                       lnormFit = TRUE,
                                                       normDeg = TRUE,
                                                       normBetw = TRUE)

summary(props_pears_SB_timepoint3_allConnections)
summary(props_pears_SB_timepoint3_allConnections)
write.table(summary(props_pears_SB_timepoint3_allConnections)[1],  file("microNetProps_summary_SB_timepoint3.txt"))
write.table(summary(props_pears_SB_timepoint3_allConnections)[4],  file("microNetProps_hubs_SB_timepoint3.txt"))

feature_nodes_periodo = props_pears_SB_timepoint3_allConnections$clustering$clust1
write.csv(feature_nodes_periodo,  file("microNetProps_clusters_SB_timepoint3.csv"))
interesting_feautre_nodes = feature_nodes_periodo[which(names(feature_nodes_periodo)%in%interesting_bacteria)]
write.csv(interesting_feautre_nodes,  file("microNetProps_interesting_clusters_SB_timepoint3.csv"))

table(feature_nodes_periodo)

#Define colors used
colors_manual_peri=c(wes_palette("GrandBudapest2")[2],"#FFFF00", wes_palette("Darjeeling2")[4], 
                     wes_palette("Darjeeling2")[2],  "#FF0000",wes_palette("Darjeeling1")[4],  
                     wes_palette("Darjeeling1")[3], "#00FF00",wes_palette("FantasticFox1")[2], 
                     "#8A2BE2", wes_palette("Moonrise2")[3], wes_palette("Moonrise3")[1], wes_palette("IsleofDogs1")[5], 
                     wes_palette("Moonrise3")[2], wes_palette("Darjeeling1")[1], wes_palette("IsleofDogs1")[1], 
                     wes_palette("Moonrise3"),wes_palette("Darjeeling1"))
#Plot the network
png("SB_timepoint3_AllConnections_MedianAbundance100Manual_labda01_nlabda30_NoPulsar.png", width = 10, height = 10, units = "in", res=600)
plot(props_pears_SB_timepoint3_allConnections, 
     nodeColor = "feature", 
     featVecCol = feature_nodes_periodo,
     colorVec = colors_manual_peri,
     nodeSize = "CSS",
     cexNodes = 2,
     nodeSizeSpread = 2, 
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     repulsion = 1.1, 
     rmSingles = "none",
     nodeFilter = "none", 
     nodeFilterPar = 1, 
     nodeTransp = 10, 
     hubTransp = 10, 
     highlightHubs = T,
     cexLabels=4,
     showTitle = TRUE,
     title1 = "",
     cexTitle = 10)
dev.off()

#Differential network 
net_spieceasi_SB_difference <- netConstruct(data = t(abundances_SB_timepoint3_4),
                                            data2 = t(abundances_SB_timepoint1_4),
                                                           dataType = "counts",
                                                           filtTax = c("none"),
                                                           measure = "spieceasi", 
                                                           measurePar=list(method = "mb", nlambda = 10, 
                                                                           lambda.min.ratio = 0.1),
                                                           jointPrepro = NULL, 
                                                           filtSamp = "none", filtSampPar = NULL,  zeroMethod = "none",
                                                           zeroPar = NULL, normMethod = "none", normPar = NULL, sparsMethod = "none",
                                                           cores = 20,  dissFunc = "signed",  weighted = TRUE, verbose = 2, 
                                                           seed = 123456)
#Analyze the network
props_pears_SB_difference <- netAnalyze(net_spieceasi_SB_difference, 
                                                       centrLCC = FALSE,
                                                       avDissIgnoreInf = TRUE,
                                                       sPathNorm = FALSE,
                                                       clustMethod = "cluster_fast_greedy",
                                                       hubPar = c("between", "degree"),
                                                       hubQuant = 0.85,
                                                       lnormFit = TRUE,
                                                       normDeg = TRUE,
                                                       normBetw = TRUE)

summary(props_pears_SB_difference)

png("SB_difference_AllConnections_MedianAbundance100Manual_labda01_nlabda10_NoPulsar.png", width = 10, height = 10, units = "in", res=600)
plot(props_pears_SB_difference, 
    sameLayout = TRUE, 
    layoutGroup = 1,
    rmSingles = "inboth", 
    nodeSize = "CSS",
    cexNodes = 2,
    nodeSizeSpread = 2, 
    edgeTranspLow = 30, 
    edgeTranspHigh = 10,
    repulsion = 1.1, 
    nodeFilter = "none", 
    nodeFilterPar = 1, 
    nodeTransp = 10, 
    hubTransp = 10, 
    highlightHubs = T,
    cexLabels=4,
    showTitle = TRUE,
    title1 = "",
    cexTitle = 0.8,
    groupNames = c("Εστιασμένη Περιοδοντίτιδα", "Γενικευμένη Περιοδοντίτιδα"),
    hubBorderCol  = "gray40")
dev.off()


write.table(summary(props_pears_SB_difference)[1],  file("microNetProps_summary_SB_difference.txt"))
write.table(summary(props_pears_SB_difference)[4],  file("microNetProps_hubs_SB_difference.txt"))


#Lets start with the SP tissue in timepoint 1 
abundances_SP_timepoint1=abundances[,which(colnames(abundances)%in%metadata_timepoint1[which(metadata_timepoint1$Sample_type_abbreviation=="SP"),"Sample_name"])]
rownames(abundances_SP_timepoint1)=abundances$Name
metadata_SP = metadata_timepoint1[which(metadata_timepoint1$Sample_type_abbreviation=="SP"),]

#Find the relative frequency of microbiomes in each group
abundances_SP_timepoint1_2= apply(abundances_SP_timepoint1,2,function(x) x/sum(x))
rownames(abundances_SP_timepoint1_2)=rownames(abundances_SP_timepoint1)
abundances_SP_timepoint1_3=abundances_SP_timepoint1_2
abundances_SP_timepoint1_3= cbind(abundances_SP_timepoint1_2,median = unlist(apply(abundances_SP_timepoint1_2, 1, median)))
abundances_SP_timepoint1_3 =  data.frame(cbind(abundances_SP_timepoint1_3,mean = unlist(apply(abundances_SP_timepoint1_2, 1, mean))))

abundances_SP_timepoint1_3[which(rownames(abundances_SP_timepoint1_3)=="Porphyromonas gingivalis"),]
abundances_SP_timepoint1_3[which(rownames(abundances_SP_timepoint1_3)=="Tannerella forsythia"),]

#Find which microbiomes we will keep
keep_SP_timepoint1_median = rownames(abundances_SP_timepoint1_3[order(abundances_SP_timepoint1_3$median, decreasing = TRUE), ])[1:100]
keep_SP_timepoint1_mean = rownames(abundances_SP_timepoint1_3[order(abundances_SP_timepoint1_3$mean, decreasing = TRUE), ])[1:100]

#Keep the 100 microbiomes with the highest median relative frequency in each group 
abundances_SP_timepoint1_4=abundances_SP_timepoint1[which(rownames(abundances_SP_timepoint1)%in%keep_SP_timepoint1_median), ]

#Create a network using Spieceasi with positive and negative connections
#samples as rows and microbiomes as columns
net_spieceasi_SP_timepoint1_allConnections <- netConstruct(t(abundances_SP_timepoint1_4),
                                                           dataType = "counts",
                                                           filtTax = c("none"),
                                                           measure = "spieceasi", 
                                                           measurePar=list(method = "mb", nlambda = 10, 
                                                                           lambda.min.ratio = 0.1),
                                                           jointPrepro = NULL, 
                                                           filtSamp = "none", filtSampPar = NULL,  zeroMethod = "none",
                                                           zeroPar = NULL, normMethod = "none", normPar = NULL, sparsMethod = "none",
                                                           cores = 20,  dissFunc = "signed",  weighted = TRUE, verbose = 2, 
                                                           seed = 123456)


#Analyze the network
props_pears_SP_timepoint1_allConnections <- netAnalyze(net_spieceasi_SP_timepoint1_allConnections, 
                                                       centrLCC = FALSE,
                                                       avDissIgnoreInf = TRUE,
                                                       sPathNorm = FALSE,
                                                       clustMethod = "cluster_fast_greedy",
                                                       hubPar = c("between", "degree"),
                                                       hubQuant = 0.85,
                                                       lnormFit = TRUE,
                                                       normDeg = TRUE,
                                                       normBetw = TRUE)

summary(props_pears_SP_timepoint1_allConnections)
summary(props_pears_SP_timepoint1_allConnections)
write.table(summary(props_pears_SP_timepoint1_allConnections)[1],  file("microNetProps_summary_SP_timepoint1.txt"))
write.table(summary(props_pears_SP_timepoint1_allConnections)[4],  file("microNetProps_hubs_SP_timepoint1.txt"))

feature_nodes_periodo = props_pears_SP_timepoint1_allConnections$clustering$clust1
write.csv(feature_nodes_periodo,  file("microNetProps_clusters_SP_timepoint1.csv"))
interesting_feautre_nodes = feature_nodes_periodo[which(names(feature_nodes_periodo)%in%interesting_bacteria)]
write.csv(interesting_feautre_nodes,  file("microNetProps_interesting_clusters_SP_timepoint1.csv"))

table(feature_nodes_periodo)

#Define colors used
colors_manual_peri=c(wes_palette("GrandBudapest2")[2],"#FFFF00", wes_palette("Darjeeling2")[4], 
                     wes_palette("Darjeeling2")[2],  "#FF0000",wes_palette("Darjeeling1")[4],  
                     wes_palette("Darjeeling1")[3], "#00FF00",wes_palette("FantasticFox1")[2], 
                     "#8A2BE2", wes_palette("Moonrise2")[3], wes_palette("Moonrise3")[1], wes_palette("IsleofDogs1")[5], 
                     wes_palette("Moonrise3")[2], wes_palette("Darjeeling1")[1], wes_palette("IsleofDogs1")[1], 
                     wes_palette("Moonrise3"),wes_palette("Darjeeling1"))
#Plot the network
png("SP_timepoint1_AllConnections_MedianAbundance100Manual_labda01_nlabda30_NoPulsar.png", width = 10, height = 10, units = "in", res=600)
plot(props_pears_SP_timepoint1_allConnections, 
     nodeColor = "feature", 
     featVecCol = feature_nodes_periodo,
     colorVec = colors_manual_peri,
     nodeSize = "CSS",
     cexNodes = 2,
     nodeSizeSpread = 2, 
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     repulsion = 1.1, 
     rmSingles = "none",
     nodeFilter = "none", 
     nodeFilterPar = 1, 
     nodeTransp = 10, 
     hubTransp = 10, 
     highlightHubs = T,
     cexLabels=4,
     showTitle = TRUE,
     title1 = "",
     cexTitle = 10)
dev.off()


#Lets continue with the SP tissue in timepoint 3
abundances_SP_timepoint3=abundances[,which(colnames(abundances)%in%metadata_timepoint3[which(metadata_timepoint3$Sample_type_abbreviation=="SP"),"Sample_name"])]
rownames(abundances_SP_timepoint3)=abundances$Name
metadata_SP = metadata_timepoint3[which(metadata_timepoint3$Sample_type_abbreviation=="SP"),]

#Find the relative frequency of microbiomes in each group
abundances_SP_timepoint3_2= apply(abundances_SP_timepoint3,2,function(x) x/sum(x))
rownames(abundances_SP_timepoint3_2)=rownames(abundances_SP_timepoint3)
abundances_SP_timepoint3_3=abundances_SP_timepoint3_2
abundances_SP_timepoint3_3= cbind(abundances_SP_timepoint3_2,median = unlist(apply(abundances_SP_timepoint3_2, 1, median)))
abundances_SP_timepoint3_3 =  data.frame(cbind(abundances_SP_timepoint3_3,mean = unlist(apply(abundances_SP_timepoint3_2, 1, mean))))

abundances_SP_timepoint3_3[which(rownames(abundances_SP_timepoint3_3)=="Porphyromonas gingivalis"),]
abundances_SP_timepoint3_3[which(rownames(abundances_SP_timepoint3_3)=="Edwardsiella ictaluri"),]

#Find which microbiomes we will keep
keep_SP_timepoint3_median = rownames(abundances_SP_timepoint3_3[order(abundances_SP_timepoint3_3$median, decreasing = TRUE), ])[1:100]
keep_SP_timepoint3_mean = rownames(abundances_SP_timepoint3_3[order(abundances_SP_timepoint3_3$mean, decreasing = TRUE), ])[1:100]

#Keep the 100 microbiomes with the highest median relative frequency in each group 
abundances_SP_timepoint3_4=abundances_SP_timepoint3[which(rownames(abundances_SP_timepoint3)%in%keep_SP_timepoint3_median), ]

#Create a network using Spieceasi with positive and negative connections
#samples as rows and microbiomes as columns
net_spieceasi_SP_timepoint3_allConnections <- netConstruct(t(abundances_SP_timepoint3_4),
                                                           dataType = "counts",
                                                           filtTax = c("none"),
                                                           measure = "spieceasi", 
                                                           measurePar=list(method = "mb", nlambda = 10, 
                                                                           lambda.min.ratio = 0.1),
                                                           jointPrepro = NULL, 
                                                           filtSamp = "none", filtSampPar = NULL,  zeroMethod = "none",
                                                           zeroPar = NULL, normMethod = "none", normPar = NULL, sparsMethod = "none",
                                                           cores = 20,  dissFunc = "signed",  weighted = TRUE, verbose = 2, 
                                                           seed = 123456)


#Analyze the network
props_pears_SP_timepoint3_allConnections <- netAnalyze(net_spieceasi_SP_timepoint3_allConnections, 
                                                       centrLCC = FALSE,
                                                       avDissIgnoreInf = TRUE,
                                                       sPathNorm = FALSE,
                                                       clustMethod = "cluster_fast_greedy",
                                                       hubPar = c("between", "degree"),
                                                       hubQuant = 0.85,
                                                       lnormFit = TRUE,
                                                       normDeg = TRUE,
                                                       normBetw = TRUE)

summary(props_pears_SP_timepoint3_allConnections)
summary(props_pears_SP_timepoint3_allConnections)
write.table(summary(props_pears_SP_timepoint3_allConnections)[1],  file("microNetProps_summary_SP_timepoint3.txt"))
write.table(summary(props_pears_SP_timepoint3_allConnections)[4],  file("microNetProps_hubs_SP_timepoint3.txt"))

feature_nodes_periodo = props_pears_SP_timepoint3_allConnections$clustering$clust1
write.csv(feature_nodes_periodo,  file("microNetProps_clusters_SP_timepoint3.csv"))
interesting_feautre_nodes = feature_nodes_periodo[which(names(feature_nodes_periodo)%in%interesting_bacteria)]
write.csv(interesting_feautre_nodes,  file("microNetProps_interesting_clusters_SP_timepoint3.csv"))

table(feature_nodes_periodo)

#Define colors used
colors_manual_peri=c(wes_palette("GrandBudapest2")[2],"#FFFF00", wes_palette("Darjeeling2")[4], 
                     wes_palette("Darjeeling2")[2],  "#FF0000",wes_palette("Darjeeling1")[4],  
                     wes_palette("Darjeeling1")[3], "#00FF00",wes_palette("FantasticFox1")[2], 
                     "#8A2BE2", wes_palette("Moonrise2")[3], wes_palette("Moonrise3")[1], wes_palette("IsleofDogs1")[5], 
                     wes_palette("Moonrise3")[2], wes_palette("Darjeeling1")[1], wes_palette("IsleofDogs1")[1], 
                     wes_palette("Moonrise3"),wes_palette("Darjeeling1"))
#Plot the network
png("SP_timepoint3_AllConnections_MedianAbundance100Manual_labda01_nlabda30_NoPulsar.png", width = 10, height = 10, units = "in", res=600)
plot(props_pears_SP_timepoint3_allConnections, 
     nodeColor = "feature", 
     featVecCol = feature_nodes_periodo,
     colorVec = colors_manual_peri,
     nodeSize = "CSS",
     cexNodes = 2,
     nodeSizeSpread = 2, 
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     repulsion = 1.1, 
     rmSingles = "none",
     nodeFilter = "none", 
     nodeFilterPar = 1, 
     nodeTransp = 10, 
     hubTransp = 10, 
     highlightHubs = T,
     cexLabels=4,
     showTitle = TRUE,
     title1 = "",
     cexTitle = 10)
dev.off()


#Differential network 
net_spieceasi_SP_difference <- netConstruct(data = t(abundances_SP_timepoint3_4),
                                            data2 = t(abundances_SP_timepoint1_4),
                                            dataType = "counts",
                                            filtTax = c("none"),
                                            measure = "spieceasi", 
                                            measurePar=list(method = "mb", nlambda = 10, 
                                                            lambda.min.ratio = 0.1),
                                            jointPrepro = NULL, 
                                            filtSamp = "none", filtSampPar = NULL,  zeroMethod = "none",
                                            zeroPar = NULL, normMethod = "none", normPar = NULL, sparsMethod = "none",
                                            cores = 20,  dissFunc = "signed",  weighted = TRUE, verbose = 2, 
                                            seed = 123456)
#Analyze the network
props_pears_SP_difference <- netAnalyze(net_spieceasi_SP_difference, 
                                        centrLCC = FALSE,
                                        avDissIgnoreInf = TRUE,
                                        sPathNorm = FALSE,
                                        clustMethod = "cluster_fast_greedy",
                                        hubPar = c("between", "degree"),
                                        hubQuant = 0.85,
                                        lnormFit = TRUE,
                                        normDeg = TRUE,
                                        normBetw = TRUE)

summary(props_pears_SP_difference)

png("SP_difference_AllConnections_MedianAbundance100Manual_labda01_nlabda10_NoPulsar.png", width = 10, height = 10, units = "in", res=600)
plot(props_pears_SP_difference, 
     sameLayout = TRUE, 
     layoutGroup = 1,
     rmSingles = "inboth", 
     nodeSize = "CSS",
     cexNodes = 2,
     nodeSizeSpread = 2, 
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     repulsion = 1.1, 
     nodeFilter = "none", 
     nodeFilterPar = 1, 
     nodeTransp = 10, 
     hubTransp = 10, 
     highlightHubs = T,
     cexLabels=4,
     showTitle = TRUE,
     title1 = "",
     cexTitle = 0.8,
     groupNames = c("Εστιασμένη Περιοδοντίτιδα", "Γενικευμένη Περιοδοντίτιδα"),
     hubBorderCol  = "gray40")
dev.off()


write.table(summary(props_pears_SP_difference)[1],  file("microNetProps_summary_SP_difference.txt"))
write.table(summary(props_pears_SP_difference)[4],  file("microNetProps_hubs_SP_difference.txt"))

