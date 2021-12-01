library(Seurat)
source("scripts/helperFunctions.R")

setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))

GE_Data <- readRDS(paste(snakemake@input[[1]], sep=""))
projectDirPath <- snakemake@params[[1]]
samples <- snakemake@params[[2]]
names <- snakemake@params[[3]]

for(i in 1:length(GE_Data)) {
  print(samples[[i]])
  print(GE_Data[[i]])
  GE_Data[[i]] <- subset(GE_Data[[i]], subset = percent.mt < samples[[i]])
  print(GE_Data[[i]])
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", names[[i]], ".mt.after.pdf", sep=""),
                   to_plot=VlnPlot(GE_Data[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
}
saveRDS(GE_Data, file=snakemake@output[[1]])