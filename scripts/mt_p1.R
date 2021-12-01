library(Seurat)
source("scripts/helperFunctions.R")

setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

names <- snakemake@params[[2]]
GE_Data <- readRDS(snakemake@input[[1]])
projectDirPath <- snakemake@params[[1]]
pattern <- snakemake@params[[3]]

for(i in 1:length(GE_Data)) {
  print(names[[i]])
  GE_Data[[i]][["percent.mt"]] <- PercentageFeatureSet(GE_Data[[i]], pattern=pattern)
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", names[[i]], ".mt.before.pdf", sep=""),
                   to_plot=VlnPlot(GE_Data[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  print(mean(GE_Data[[i]]@meta.data$percent.mt))
  print(median(GE_Data[[i]]@meta.data$percent.mt))
}
saveRDS(GE_Data, file=snakemake@output[[1]])