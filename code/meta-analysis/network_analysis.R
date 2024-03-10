# Install NetCoMi
 # devtools::install_github("stefpeschel/NetCoMi", 
 #                          dependencies = c("Depends", "Imports", "LinkingTo"),
 #                          repos = c("https://cloud.r-project.org/",
 #                                    BiocManager::repositories()))

library(NetCoMi)
library(MetBrewer)
library(data.table)
library(wesanderson)
library(dplyr)

#Read the file with abundances of samples before normalisation
data=read.delim(file="./species_nozero.csv", header=T)
data2=data.frame(t(data))
colnames(data2)=data2[1,]
data2=data2[2:5177,]

#Read the metadata file 
metadata=read.delim(file='./RNA_metadata_ALL.tab')


#Core species 
core_species=read.delim(file="All_species_core.txt")
diff_species=read.delim(file="differential_expressed_species.txt",sep=",", header=F)

#kingdom of microbiomes 
kingdom=read.delim(file="kingdom_species_armen.txt", sep=",", header=F)
bacteria=kingdom$V2[which(kingdom$V1=="Bacteria")]
bacteria = sub(" ",".",bacteria)

#Keep only the microbiomes in bacteria kingdom
data3=data2[which(rownames(data2)%in%bacteria),]

#Separate the dataframe in healthy and periodontitis and make all the values numeric
data4_periodontitis=data3[,which(metadata$Condition=="Periodontitis")]
data4_healthy=data3[,which(metadata$Condition=="Healthy")]
data4_periodontitis_2=sapply(data4_periodontitis, function(x) as.numeric(as.character(x)))
data4_healthy_2=sapply(data4_healthy, function(x) as.numeric(as.character(x)))

data4_periodontitis[1,1]
data4_periodontitis_2[1,1]
rownames(data4_periodontitis_2)=rownames(data4_periodontitis)
rownames(data4_healthy_2)=rownames(data4_healthy)

#check to see if treponema or tannerella have higher mean
rowSums(data4_periodontitis_2)[which(rownames(data4_periodontitis_2)%in%c("Treponema.denticola"))]
rowSums(data4_periodontitis_2)[which(rownames(data4_periodontitis_2)%in%c("Tannerella.forsythia"))]
summary(data4_periodontitis_2[which(rownames(data4_periodontitis_2)%in%c("Tannerella.forsythia")), ])
summary(data4_periodontitis_2[which(rownames(data4_periodontitis_2)%in%c("Treponema.denticola")), ])

#Find the relative frequency of microbiomes in each group
sum(data4_healthy_2[,1])
data4_healthy_3= apply(data4_healthy_2,2,function(x) x/sum(x))
data4_periodontitis_3= apply(data4_periodontitis_2,2,function(x){x/sum(x)})
rownames(data4_periodontitis_3)=rownames(data4_periodontitis_2)
rownames(data4_healthy_3)=rownames(data4_healthy_2)
data4_healthy_4=data4_healthy_3
data4_healthy_4= cbind(data4_healthy_3,median = unlist(apply(data4_healthy_3, 1, median)))
data4_healthy_4 =  data.frame(cbind(data4_healthy_4,mean = unlist(apply(data4_healthy_3, 1, mean))))
data4_periodontitis_4=data4_periodontitis_3
data4_periodontitis_4= cbind(data4_periodontitis_3,median = unlist(apply(data4_periodontitis_3, 1, median)))
data4_periodontitis_4 = data.frame(cbind(data4_periodontitis_4,mean = unlist(apply(data4_periodontitis_3, 1, mean))))

data4_periodontitis_4[which(rownames(data4_periodontitis_4)%in%("Tanne")), ]
#Find which microbiomes we will keep
keep_healthy_median = rownames(data4_healthy_4[order(data4_healthy_4$median, decreasing = TRUE), ])[1:100]
keep_healthy_mean = rownames(data4_healthy_4[order(data4_healthy_4$mean, decreasing = TRUE), ])[1:100]
keep_periodontitis_median = rownames(data4_periodontitis_4[order(data4_periodontitis_4$median, decreasing = TRUE), ])[1:100]
keep_periodontitis_mean = rownames(data4_periodontitis_4[order(data4_periodontitis_4$mean, decreasing = TRUE), ])[1:100]

#Keep the 100 microbiomes with the highest median relative frequency in each group 
data4_healthy_5=data4_healthy_2[which(rownames(data4_healthy_2)%in%keep_healthy_median), ]
length(which(keep_healthy_mean%in%core_species$Row.names))
length(which(keep_healthy_median%in%core_species$Row.names))

data4_periodontitis_5=data4_periodontitis_2[which(rownames(data4_periodontitis_2)%in%keep_periodontitis_median), ]
length(which(keep_periodontitis_median%in%core_species$Row.names))
length(which(keep_periodontitis_mean%in%core_species$Row.names))

#Create a network using Spieceasi with positive and negative connections
#samples as rows and microbiomes as columns
net_spieceasi_periodontitis_allConnections = netConstruct(t(data4_periodontitis_5),
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
props_pears_periodo_allConnections = netAnalyze(net_spieceasi_periodontitis_allConnections, 
                          centrLCC = FALSE,
                          avDissIgnoreInf = TRUE,
                          sPathNorm = FALSE,
                          clustMethod = "cluster_fast_greedy",
                          hubPar = c("between", "degree"),
                          hubQuant = 0.85,
                          lnormFit = TRUE,
                          normDeg = TRUE,
                          normBetw = TRUE)

summary(props_pears_periodo_allConnections)


feature_nodes_periodo = props_pears_periodo_allConnections$clustering$clust1
table(feature_nodes_periodo)

#Define colors used
colors_manual_peri=c(wes_palette("GrandBudapest2")[2],"#FFFF00", wes_palette("Darjeeling2")[4], 
                     wes_palette("Darjeeling2")[2],  "#FF0000",wes_palette("Darjeeling1")[4],  
                     wes_palette("Darjeeling1")[3], "#00FF00",wes_palette("FantasticFox1")[2], 
                     "#8A2BE2", wes_palette("Moonrise2")[3], wes_palette("Moonrise3")[1], wes_palette("IsleofDogs1")[5], 
                     wes_palette("Moonrise3")[2], wes_palette("Darjeeling1")[4], wes_palette("IsleofDogs1")[1])
#Plot the network
png("periodontitis_AllConnections_MedianAbundance100Manual_labda01_nlabda10_NoPulsar.png", width = 20, height = 20, units = "in", res=600)
plot(props_pears_periodo_allConnections, 
     nodeColor = "feature", 
     featVecCol = feature_nodes_periodo,
     colorVec = colors_manual_peri,
     nodeSize = "CSS",
     cexNodes = 2,
     nodeSizeSpread = 2, 
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     repulsion = 1.2, 
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

#Network for healthy
net_spieceasi_healthy_allConnections = netConstruct(t(data4_healthy_5),
                                                           dataType = "counts",
                                                           filtTax = c("none"),
                                                           measure = "spieceasi",
                                                           measurePar=list(method = "mb", nlambda = 10, 
                                                                     lambda.min.ratio = 0.1),
                                                           jointPrepro = NULL, 
                                                           filtSamp = "none", filtSampPar = NULL,  zeroMethod = "none",
                                                           zeroPar = NULL, normMethod = "none", normPar = NULL, sparsMethod = "none",
                                                           cores = 20, dissFunc = "signed",  weighted = TRUE, verbose = 2, 
                                                           seed = 123456)

props_pears_healthy_allConnections = netAnalyze(net_spieceasi_healthy_allConnections, 
                                  centrLCC = FALSE,
                                  avDissIgnoreInf = TRUE,
                                  sPathNorm = FALSE,
                                  clustMethod = "cluster_fast_greedy",
                                  hubPar = c("between", "degree"),
                                  hubQuant = 0.85,
                                  lnormFit = TRUE,
                                  normDeg = TRUE,
                                  normBetw = TRUE)

# summary(props_pears_healthy_positiveConnections)
summary(props_pears_healthy_allConnections)

summary(props_pears_periodo_allConnections)

feature_nodes_healthy = props_pears_healthy_allConnections$clustering$clust1
table(feature_nodes_healthy)


colors_manual_healthy=c(wes_palette("Moonrise2")[1], wes_palette("Darjeeling2")[1], "#FFFF00",
                        wes_palette("Darjeeling2")[4], "#FF0000" , "#00FF00", 
                        wes_palette("Darjeeling1")[4],wes_palette("Darjeeling2")[2], wes_palette("Moonrise2")[3], 
                        wes_palette("Moonrise1")[3], wes_palette("Moonrise3")[1], "#8A2BE2",
                        wes_palette("IsleofDogs1")[5])
#Plot the network
png("healthy_AllConnections_Mean100AbundantManual_labda01_nlabda10_NoPulsar.png", width = 20, height = 20, units = "in", res=600)
plot(props_pears_healthy_allConnections, 
     nodeColor = "feature", 
     featVecCol = feature_nodes_healthy,
     colorVec = colors_manual_healthy,
     nodeSize = "CSS",
     cexNodes = 2,
     nodeSizeSpread = 2, 
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     repulsion = 1.2, 
     rmSingles = "none",
     nodeFilter = "none", 
     nodeFilterPar = 1, 
     nodeTransp = 10, 
     hubTransp = 10, 
     highlightHubs = T,
     cexLabels=4,
     showTitle = TRUE,
     title1 = "",
     cexTitle=10)
dev.off()

