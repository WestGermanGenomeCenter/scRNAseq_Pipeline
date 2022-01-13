library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

GE.Data <- readRDS(snakemake@input[[1]])
project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
condition <- snakemake@params[[3]]
metaIdents <- unique(snakemake@params[[4]])
txtContent <- c("finishedDGE...")

Idents(GE.Data) <- paste("integrated_snn_res", snakemake@params[[5]], sep=".")
len <- length(levels(GE.Data))-1 #clusters named from 0 to n-1
for(i in 0:len) {
  #Subset Seurat Object to cluster to be analyzed
  cluster <- subset(GE.Data, idents=i)
  #Set identity to condition
  cluster <- SetIdent(cluster, value=condition)
  #Find Markers between both conditions
  if(length(levels(cluster)) == length(metaIdents)) {
    cluster.conditions <- SplitObject(cluster, split.by=condition)
    if(length(colnames(cluster.conditions[[1]])) < 3) {
      txtContent <- paste(txtContent, paste("Cluster", i, "was skipped since", levels(cluster.conditions[[1]])[[1]], "has less than 3 cells", sep=" "), sep="\n")
    } else if(length(colnames(cluster.conditions[[2]])) < 3) {
      txtContent <- paste(txtContent, paste("Cluster", i, "was skipped since", levels(cluster.conditions[[2]])[[1]], "has less than 3 cells", sep=" "), sep="\n")
    } else {
      cluster.markers <- FindMarkers(cluster, only.pos=F, min.pct=0.25, logfc.threshold=0.25, ident.1=metaIdents[[1]], ident.2=metaIdents[[2]])
      #plot_in_terminal(plotname=paste("../projects/", projectDirName, "plots/", i, ".cluster.dge.pdf", sep=""),
      #               to_plot=DoHeatmap(cluster, features=cluster.markers))
      write.csv(cluster.markers, file=paste(projectDirPath, "csv/", project, ".compareGE_inCluster", i, ".markers.csv", sep=""))
    }
  } else {
    txtContent <- paste(txtContent, paste("Cluster", i, "was skipped since it was made up of only", levels(cluster)[[1]], sep=" "), sep="\n")
  }
  rm(cluster)
}


connection <- file(paste(projectDirPath, "csv/finishedDGE.txt", sep=""))
writeLines(txtContent, connection)
close(connection)
