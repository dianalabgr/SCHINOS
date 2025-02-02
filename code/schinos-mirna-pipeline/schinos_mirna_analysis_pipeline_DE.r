## Required packages installation
packages <- c("parallel","optparse","RColorBrewer")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos='http://cran.us.r-project.org')  
}

## CMD Argument Input
library("optparse", quietly = TRUE, warn.conflicts = FALSE)

option_list <- list(
  make_option(c("-d", "--difexpression"),
              help="Provide the path to the condition table for the \n\t\tdifferential expression analysis.\n\t\tIt should contain two columns with their headers,\n\t\tthe first being the sample column containing all \n\t\tthe sample names and the second should be the \n\t\tcondition column, containing a condition for each sample. \n\t\tNOTICE: The different conditions must be exactly two in number.",
              metavar="Path")
)
pargv = parse_args(OptionParser(usage = "\n%prog [options] input_directory output_directory number_of_cores", option_list = option_list, description = "\nApply Differential Expression analysis to the analysis results."), positional_arguments = 3)

input_directory = pargv$args[[1]]
output_directory = pargv$args[[2]]
core_num = pargv$args[[3]]

if (substr(input_directory,nchar(input_directory),nchar(input_directory)) == "/") {
  input_directory = substr(input_directory, 1, nchar(input_directory)-1)
}

library("parallel", quietly = TRUE, warn.conflicts = FALSE)

## Differential Expression Analysis
if (exists("difexpression")) {
  cat("DE ANALYSIS START TIME: ")
  system("date")
  cat("=> INITIATING DIFFERENTIAL EXPRESSION ANALYSIS... \n")
  bioconductor_packages <- c("DESeq2","pcaExplorer")
  if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(bioconductor_packages)
  }
  suppressMessages(library("DESeq2", quietly = TRUE, warn.conflicts = FALSE))
  suppressMessages(library("pcaExplorer", quietly = TRUE, warn.conflicts = FALSE))
  
  dir.create(output_directory)
  setwd(output_directory)
  
  ## Read-in conditions table
  colData = read.csv(file = difexpression, stringsAsFactors = F)
  colData[[colnames(colData)[2]]] = factor(colData[[colnames(colData)[2]]])
  if (nlevels(colData[[colnames(colData)[2]]]) != 2) {
    stop(paste0("The conditions provided for the Differential Expression of the samples are not 2.\nPlease change the input table accordingly in order to provide exactly 2 different conditions.\n\nConditions given: ",toString(levels(factor(colData[[colnames(colData)[2]]]))),"\n"))
  }
  
  ## Create a counts table
  for(sample in colData[,1]) {
    sample_counts_file = Sys.glob(file.path(input_directory,paste0(sample,"*"),paste0(sample,"*_Counts.txt")))
    if (length(sample_counts_file) > 0 && file.exists(sample_counts_file)) {
      sample_counts = read.table(sample_counts_file, sep="\t", col.names = c("miRNA_ID",paste0(sample,"_raw_counts"),as.character(sample),paste0(sample,"_log2")), header = FALSE, skip = 1 )
      
      if (!exists("master_table")){
        master_table = sample_counts[,c(1,2)]
      }else {
        master_table = merge(master_table,sample_counts[,c(1,2)],by = "miRNA_ID")
      }
    }
  }
  
  write.csv(master_table, file = "raw_counts_table.csv",row.names=FALSE)
  
  ## Read-in the counts table
  countFilePath = "raw_counts_table.csv"
  countData = read.table(file = countFilePath, header = TRUE, sep = ",", row.names = 1)
  countData = round(countData,0)
  colnames(countData) = sub("_raw_counts","",colnames(countData))
  
  group_name = colnames(colData)[2]
  group1 = levels(factor(colData[[colnames(colData)[2]]]))[1]
  group2 = levels(factor(colData[[colnames(colData)[2]]]))[2]
  
  ## Import the counts table and run DESeq2 for the DE Analysis
  dataset = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = as.formula(paste0("~",group_name)))
  
  dds = DESeq(dataset)
  
  result = results(dds, contrast=c(group_name,group1,group2))
  result = result[complete.cases(result),] # Remove any rows with NA
  
  ## Write-out the results
  write.table(counts(dds, normalized=TRUE), file="DE_normalized_counts.txt", sep="\t", quote=F, col.names=NA)
  write.csv(as.data.frame(result),file = "DE_full_results.csv")
  writeLines(c("\t\t\t\tDIFFERENTIAL EXPRESSION REPORT\n\n","\t\t\tSUMMARY:"),"DE_report.txt")
  capture.output(summary(result),file = "DE_report.txt", append = TRUE)
  
  ## DE Plots
  # PCA
  cat(">> PLOTTING PCA... \n")
  rld = rlogTransformation(dds, blind = TRUE)
  pcaobj = prcomp(t(assay(rld)))
  pdf("DE_PCA.pdf")
  print(pcaplot(rld, intgroup = group_name,title = paste0("Principal Component Analysis on ",group_name)))
  print(pcascree(pcaobj,type = "pev",title = "Proportion of Explained Variance per Principal Component"))
  dev.off()
  
  # MA
  cat(">> PLOTTING MA-PLOT... \n")
  png(filename="DE_MAplot.png")
  plotMA(result, main=paste0("MA-plot on ",group_name,": ",group1," vs. ",group2), ylim=c(-5,5))
  dev.off()
  
  ## Extract the top n up-regulated and the top n down-regulated miRNAs, sorted by adjusted p-value
  cat(">> EXTRACTING TOP RESULTS... \n")
  n = 100
  resOrdered = result[order(result$padj),]
  topResults = rbind(resOrdered[resOrdered[,"log2FoldChange"] > 0, ][1:n,], resOrdered[resOrdered[,"log2FoldChange"] < 0, ][n:1,])
  write.csv(as.data.frame(topResults),file = paste0("DE_top",n,"_by_padj_results.csv"))

  ## Report top and bottom 5 miRNAs
  de_report = readLines("DE_report.txt")
  writeLines(c(de_report," ","\t\t\tTOP AND BOTTOM 5 miRNAs BASED ON ADJUSTED P-VALUE:"),"DE_report.txt")
  capture.output(topResults[c(1:5,(2*n-4):(2*n)), c("baseMean","log2FoldChange","padj")],file = "DE_report.txt", append = TRUE)

  ## Plot count comparisons per top miRNA
  cat(">> PLOTTING COUNT COMPARISON PLOTS PER miRNA... \n")
  pdf(paste0("top",n,"_miRNA_count_comparison_plots.pdf"))
  for (gn in row.names(topResults)) {
    plotCounts(dds, gene = gn, intgroup = group_name, pch = 19)
  }
  dev.off()
  
  ## DE Plot Heatmaps
  # Plot Top 100 miRNAs heatmap and cluster by similarity
  cat(">> PLOTTING HEATMAPS... \n")
  suppressMessages(library("RColorBrewer", quietly = TRUE, warn.conflicts = FALSE))
  hmcol = brewer.pal(11,"RdBu")
  nCounts = counts(dds, normalized = TRUE)
  png(filename=paste0("DE_top",n,"_by_padj_heatmap.png"), width = 8, height = 8, units = 'in', res = 600)
  heatmap(as.matrix(nCounts[row.names(topResults),]), Rowv = NA, col = hmcol, mar = c(10,2))
  dev.off()
  
  # Plot Top and Bottom 25 miRNAs heatmap and cluster by similarity
  m = 25
  png(filename=paste0("DE_top",m,"_by_padj_heatmap.png"), width = 6, height = 6, units = 'in', res = 600)
  heatmap(as.matrix(nCounts[row.names(topResults)[c(1:m,(n-m+1):n)], ]), Rowv = NA, col = hmcol, mar = c(10,2))
  dev.off()
  
  cat("=> DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED \n")
  cat("DE ANALYSIS END TIME: ")
  system("date")
}
