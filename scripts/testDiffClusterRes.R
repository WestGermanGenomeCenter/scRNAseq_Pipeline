library(Seurat)
library(future)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

project.name <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
dimens <- snakemake@params[[3]]
resolutions <- snakemake@params[[4]]
cores <- snakemake@params[[5]]
GE.integrated <- readRDS(snakemake@input[[1]])
GE.integrated <- FindNeighbors(GE.integrated, dims=1:dimens)

plan("multiprocess", workers=cores)
GE.integrated <- FindClusters(GE.integrated, resolution=resolutions)
for(i in 1:length(resolutions)) {
  resolution <- resolutions[[i]]
  #Future when multiples are use
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", project.name, ".res_", resolution, ".clusteredDimPlot.pdf", sep=""),
                  to_plot=DimPlot(GE.integrated, group.by=paste("integrated_snn_res", resolution, sep=".")))
}
