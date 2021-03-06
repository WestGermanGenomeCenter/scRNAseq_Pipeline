library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

project.name <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
dimens <- snakemake@params[[3]]
resolution <- snakemake@params[[4]]
condition <- snakemake@params[[5]]
GE.integrated <- readRDS(snakemake@input[[1]])
GE.integrated <- FindNeighbors(GE.integrated, dims=1:dimens)

GE.integrated <- FindClusters(GE.integrated, resolution=resolution)
plot_in_terminal(plotname=paste(projectDirPath, "plots/", project.name, ".labeledDimPlot.pdf", sep=""),
                 to_plot=DimPlot(GE.integrated, label=T))
print(condition)
for(i in 1:length(condition)) {
    print(condition[[i]])
    plot_in_terminal(plotname=paste(projectDirPath, "plots/", project.name, ".", condition[[i]], "DimPlot.pdf", sep=""),
                 to_plot=DimPlot(GE.integrated,split.by=condition[[i]]))
}

saveRDS(GE.integrated, file=snakemake@output[[1]])
