library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))
GE_Data <- readRDS(snakemake@input[[1]])
names <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]

#Doublet Removal
for(i in 1:length(GE_Data)) {
  GE_Data[[i]] <- SCTransform(GE_Data[[i]], vars.to.regress="percent.mt", verbose=FALSE)
  GE_Data[[i]] <- RunPCA(GE_Data[[i]], verbose=FALSE)
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", names[[i]], ".ElbowPlot.pdf", sep=""),
                   to_plot=ElbowPlot(GE_Data[[i]]))
}
saveRDS(GE_Data, file=snakemake@output[[1]])