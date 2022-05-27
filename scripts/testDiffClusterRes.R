library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

project.name <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
dimens <- snakemake@params[[3]]
resolutions <- snakemake@params[[4]]
GE.integrated <- readRDS(snakemake@input[[1]])
GE.integrated <- FindNeighbors(GE.integrated, dims=1:dimens)

for(i in 1:length(resolutions)) {
  resolution <- resolutions[[i]]
  #Future when multiples are use
  GE.integrated <- FindClusters(GE.integrated, resolution=resolution)
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", project.name, ".res_", resolution, ".clusteredDimPlot.pdf", sep=""),
                  to_plot=DimPlot(GE.integrated))
}