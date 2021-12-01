library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))

project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
dimens <- snakemake@params[[3]]
GE.215.integrated <- readRDS(snakemake@input[[1]])

GE.215.integrated <- RunUMAP(GE.215.integrated, dims=1:dimens)#:30)
plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".umappedDimPlot.pdf", sep=""),
                 to_plot=DimPlot(GE.215.integrated, group.by="sample"))
plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".umappedElbowPlot.pdf", sep=""),
                 to_plot=ElbowPlot(GE.215.integrated, ndims=dimens))

saveRDS(GE.215.integrated, file=snakemake@output[[1]])
