library(Seurat)
source("scripts/helperFunctions.R")

setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

sample.name <- snakemake@params[[2]]
GE.data <- readRDS(snakemake@input[[1]])
projectDirPath <- snakemake@params[[1]]
pattern <- snakemake@params[[3]]

print(sample.name)
GE.data[["percent.mt"]] <- PercentageFeatureSet(GE.data, pattern=pattern)
plot_in_terminal(plotname=paste(projectDirPath, "plots/", sample.name, ".mt.before.pdf", sep=""),
                 to_plot=VlnPlot(GE.data, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3))
print(mean(GE.data@meta.data$percent.mt))
print(median(GE.data@meta.data$percent.mt))
saveRDS(GE.data, file=snakemake@output[[1]])
