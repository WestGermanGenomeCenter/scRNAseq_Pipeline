library("Seurat")
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

GE <- readRDS(snakemake@input[[1]])
project.name <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
assay.name <- snakemake@params[[3]]

features <- rownames(GE@assays[[assay.name]]@data)

DefaultAssay(GE) <- assay.name

for(i in 1:length(features)) {
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", features[[i]], ".adtFeatureProtein.pdf", sep=""),
                   to_plot=FeaturePlot(GE, features=features[[i]]))
}

connection <- file(paste(projectDirPath, "plots/finishedFeatures.txt", sep=""))
writeLines(c("finished"), connection)
close(connection)
