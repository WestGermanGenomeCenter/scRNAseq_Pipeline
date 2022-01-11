library(Seurat)
source("scripts/helperFunctions.R")

setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))

GE.data <- readRDS(paste(snakemake@input[[1]], sep=""))
projectDirPath <- snakemake@params[[1]]
samples <- snakemake@params[[2]]
sample.names <- snakemake@params[[3]]

for(i in 1:length(GE.data)) {
  print(samples[[i]])
  print(GE.data[[i]])
  GE.data[[i]] <- subset(GE.data[[i]], subset= percent.mt<samples[[i]])
  print(GE.data[[i]])
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", sample.names[[i]], ".mt.after.pdf", sep=""),
                   to_plot=VlnPlot(GE.data[[i]], features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3))
}
saveRDS(GE.data, file=snakemake@output[[1]])