library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
dimens <- snakemake@params[[3]]
resolution <- snakemake@params[[4]]
condition <- snakemake@params[[5]]
GE.215.integrated <- readRDS(snakemake@input[[1]])
GE.215.integrated <- FindNeighbors(GE.215.integrated, dims = 1:dimens)

GE.215.integrated <- FindClusters(GE.215.integrated, resolution=resolution)
plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".labeledDimPlot.pdf", sep=""),
                 to_plot=DimPlot(GE.215.integrated, label=T))
plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".timedDimPlot.pdf", sep=""),
                 to_plot=DimPlot(GE.215.integrated,split.by=condition))

saveRDS(GE.215.integrated, file=snakemake@output[[1]])
