library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))

project.name <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
dimens <- snakemake@params[[3]]
GE <- readRDS(snakemake@input[[1]])

GE <- RunUMAP(GE, dims=1:dimens)#:30)
plot_in_terminal(plotname=paste(projectDirPath, "plots/", project.name, ".umappedDimPlot.pdf", sep=""),
                 to_plot=DimPlot(GE, group.by="sample"))
plot_in_terminal(plotname=paste(projectDirPath, "plots/", project.name, ".umappedElbowPlot.pdf", sep=""),
                 to_plot=ElbowPlot(GE, ndims=dimens))

saveRDS(GE, file=snakemake@output[[1]])
