library(Seurat)
source("scripts/helperFunctions.R")

setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

sample.names <- snakemake@params[[2]]
GE.data <- readRDS(snakemake@input[[1]])
projectDirPath <- snakemake@params[[1]]
pattern <- snakemake@params[[3]]

for(i in 1:length(GE.data)) {
  print(sample.names[[i]])
  GE.data[[i]][["percent.mt"]] <- PercentageFeatureSet(GE.data[[i]], pattern=pattern)
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", sample.names[[i]], ".mt.before.pdf", sep=""),
                   to_plot=VlnPlot(GE.data[[i]], features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3))
  print(mean(GE.data[[i]]@meta.data$percent.mt))
  print(median(GE.data[[i]]@meta.data$percent.mt))
}
saveRDS(GE.data, file=snakemake@output[[1]])