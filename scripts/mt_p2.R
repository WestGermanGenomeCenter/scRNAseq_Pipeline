library(Seurat)
source("scripts/helperFunctions.R")

setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))

GE.data <- readRDS(paste(snakemake@input[[1]], sep=""))
projectDirPath <- snakemake@params[[1]]
sample.name <- snakemake@params[[2]]
sample <- list(snakemake@params[[3]]) #yes this is necessary, otherwise the subset won't work

print(sample)
print(GE.data)
GE.data <- subset(GE.data, subset= percent.mt<sample[[1]])
print(GE.data)
plot_in_terminal(plotname=paste(projectDirPath, "plots/", sample.name, ".mt.after.pdf", sep=""),
                 to_plot=VlnPlot(GE.data, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3))
saveRDS(GE.data, file=snakemake@output[[1]])
