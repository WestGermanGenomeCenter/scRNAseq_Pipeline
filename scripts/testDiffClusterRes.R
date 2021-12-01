resolutions <- snakemake@params[[4]]

library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
dimens <- snakemake@params[[3]]
GE.215.integrated <- readRDS(snakemake@input[[1]])
GE.215.integrated <- FindNeighbors(GE.215.integrated, dims = 1:dimens)

for(i in 1:length(resolutions)) {
  resolution <- resolutions[[i]]
  GE.215.integrated <- FindClusters(GE.215.integrated, resolution=resolution)
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".res_", resolution, ".clusteredDimPlot.pdf", sep=""),
                  to_plot=DimPlot(GE.215.integrated))
}

connection <- file(paste(projectDirPath, "plots/finished.txt", sep=""))
writeLines(c("finished"), connection)
close(connection)