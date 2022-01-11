library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))
GE.data <- readRDS(snakemake@input[[1]])
sample.names <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]

#Doublet Removal
for(i in 1:length(GE.data)) {
  GE.data[[i]] <- SCTransform(GE.data[[i]], vars.to.regress="percent.mt", verbose=FALSE)
  GE.data[[i]] <- RunPCA(GE.data[[i]], verbose=FALSE)
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", sample.names[[i]], ".ElbowPlot.pdf", sep=""),
                   to_plot=ElbowPlot(GE.data[[i]]))
}
saveRDS(GE.data, file=snakemake@output[[1]])