library("Seurat")
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

GE.215.integrated <- readRDS(snakemake@input[[1]])
project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
assayname <- snakemake@params[[3]]

features <- rownames(GE.215.integrated@assays[[assayname]]@data)

DefaultAssay(GE.215.integrated) <- assayname

for(i in 1:length(features)) {
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", features[[i]], ".adtFeatureProtein.pdf", sep=""),
                   to_plot=FeaturePlot(GE.215.integrated, features=features[[i]]))
  #gene <- str_remove(features[[i]], "adt_TotalSeqB-")
  #plot_in_terminal(plotname=paste("../projects/", projectDirName, "plots/", gene, ".adtFeatureRNA.pdf", sep=""),
  #                 to_plot=FeaturePlot(GE.215.integrated, features=features[[i]]))
}

connection <- file(paste(projectDirPath, "plots/finishedFeatures.txt", sep=""))
writeLines(c("finished"), connection)
close(connection)
