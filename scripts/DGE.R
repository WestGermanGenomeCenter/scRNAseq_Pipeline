library(Seurat)
#library(future)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

GE.Data <- readRDS(snakemake@input[[1]])
project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
condition <- snakemake@params[[3]]
#cores <- snakemake@params[[5]]
print(condition)
#metaIdents <- unlist(unique(GE.Data@meta.data[[condition]]))
#print(metaIdents)
txtContent <- c("finishedDGE...")

Idents(GE.Data) <- paste("integrated_snn_res", snakemake@params[[4]], sep=".")
len <- length(levels(GE.Data))-1 #clusters named from 0 to n-1
for(i in 0:len) {
  #Subset Seurat Object to cluster to be analyzed
  cluster <- subset(GE.Data, idents=i)
  #Set identity to condition
  cluster <- SetIdent(cluster, value=condition)
  #Find Markers between both conditions
  clusterMetaIdents <- levels(cluster)
  print(length(clusterMetaIdents))
  print(clusterMetaIdents)
  if(length(clusterMetaIdents) > 1) { #iff cluster has more than one condition
    condition.combi <- combn(clusterMetaIdents, 2)
    print(condition.combi)
    print(ncol(condition.combi))
    for(j in 1:ncol(condition.combi)) {
      condition1 <- condition.combi[1, j]
      condition2 <- condition.combi[2, j]
      print(paste(condition1, condition2))
      print(paste(length(colnames(subset(cluster, idents=condition1))), ncol(subset(cluster, idents=condition2))))
      print(paste(length(colnames(subset(cluster, idents=condition2))), ncol(subset(cluster, idents=condition2))))
      if(length(colnames(subset(cluster, idents=condition1))) < 3) {
        txtContent <- paste(txtContent, paste("Cluster", i, "was skipped since", condition1, "has less than 3 cells", sep=" "), sep="\n")
      } else if(length(colnames(subset(cluster, idents=condition2))) < 3) {
        txtContent <- paste(txtContent, paste("Cluster", i, "was skipped since", condition2, "has less than 3 cells", sep=" "), sep="\n")
      } else {
        #plan("multiprocess", workers=cores)
        cluster.markers <- FindMarkers(cluster, only.pos=F, min.pct=0.25, logfc.threshold=0.25, ident.1=condition1, ident.2=condition2, verbose=FALSE)
        #plot_in_terminal(plotname=paste("../projects/", projectDirName, "plots/", i, ".cluster.", condition1, ".", condition2, ".dge.pdf", sep=""),
        #                 to_plot=DoHeatmap(cluster, features=cluster.markers))
        write.csv(cluster.markers, file=paste(projectDirPath, "csv/", condition, "/", project, ".compare.", condition1, ".", condition2, ".cluster", i, ".markers.csv", sep=""))
      }
    }
  } else {
      txtContent <- paste(txtContent, paste("Cluster", i, "was skipped since it was made up of only", levels(cluster)[[1]], sep=" "), sep="\n")
  }
  rm(cluster)
}


connection <- file(snakemake@output[[1]])
writeLines(txtContent, connection)
close(connection)
