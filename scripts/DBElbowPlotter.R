library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))
GE.data <- readRDS(snakemake@input[[1]])
sample.name <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]

#Doublet Removal
GE.data <- SCTransform(GE.data, vars.to.regress="percent.mt", verbose=FALSE)
GE.data <- RunPCA(GE.data, verbose=FALSE)
plot_in_terminal(plotname=paste(projectDirPath, "plots/", sample.name, ".ElbowPlot.pdf", sep=""),
                 to_plot=ElbowPlot(GE.data))
saveRDS(GE.data, file=snakemake@output[[1]])