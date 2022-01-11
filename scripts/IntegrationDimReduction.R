library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

project.name <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
maxRAM <- snakemake@params[[3]]

GE.data <- readRDS(snakemake@input[[1]])
print(Idents(GE.data[[1]]))
print(length(Idents(GE.data[[1]])))
if(length(GE.data) > 1) {
  GE.features <- SelectIntegrationFeatures(GE.data, nfeatures=3000)
  options(future.globals.maxSize=maxRAM)
  GE.data <- PrepSCTIntegration(object.list=GE.data, anchor.features=GE.features, verbose=FALSE)
  GE.anchors <- FindIntegrationAnchors(object.list=GE.data, normalization.method="SCT", anchor.features=GE.features, verbose=FALSE)
  GE.integrated <- IntegrateData(anchorset=GE.anchors, normalization.method="SCT", verbose=FALSE)
  GE.integrated <- RunPCA(GE.integrated, verbose=FALSE)
} else {
  GE.integrated <- GE.data[[1]]
  GE.integrated <- RunPCA(GE.integrated, verbose=FALSE)
}
plot_in_terminal(plotname=paste(projectDirPath, "plots/", project.name, ".integratedElbowPlot.pdf", sep=""),
                 to_plot=ElbowPlot(GE.integrated))
saveRDS(GE.integrated, file=snakemake@output[[1]])
