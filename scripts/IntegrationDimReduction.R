library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
maxRAM <- snakemake@params[[3]]

GE.215.list <- readRDS(snakemake@input[[1]])
print(Idents(GE.215.list[[1]]))
print(length(Idents(GE.215.list[[1]])))
if(length(GE.215.list) > 1) {
  GE.215.features <- SelectIntegrationFeatures(GE.215.list, nfeatures=3000)
  options(future.globals.maxSize=maxRAM)
  GE.215.list <- PrepSCTIntegration(object.list=GE.215.list, anchor.features=GE.215.features, verbose=FALSE)
  GE.215.anchors <- FindIntegrationAnchors(object.list=GE.215.list, normalization.method="SCT", anchor.features=GE.215.features, verbose = FALSE)
  GE.215.integrated <- IntegrateData(anchorset=GE.215.anchors, normalization.method="SCT", verbose=FALSE)
  GE.215.integrated <- RunPCA(GE.215.integrated, verbose=FALSE)
} else {
  GE.215.integrated <- GE.215.list[[1]]
  GE.215.integrated <- RunPCA(GE.215.integrated, verbose=FALSE)
}
plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".integratedElbowPlot.pdf", sep=""),
                 to_plot=ElbowPlot(GE.215.integrated))
saveRDS(GE.215.integrated, file=snakemake@output[[1]])
