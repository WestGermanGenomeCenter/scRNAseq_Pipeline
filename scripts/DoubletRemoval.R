library(Seurat)
library(DoubletFinder)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))
GE.data <- readRDS(snakemake@input[[1]])
sample.names <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
expectedDoubletPercent <- snakemake@params[[3]]
dims_PCs <- snakemake@params[[4]]
pN <- 0.25
  
for(i in 1:length(GE.data)) {
  GE.data[[i]] <- RunUMAP(GE.data[[i]], dims=1:dims_PCs[[i]], verbose=FALSE)
  doublet_paramSweep <- paramSweep_v3(GE.data[[i]], PCs=1:dims_PCs[[i]], sct=TRUE)
  doublet_sweep.stats <- summarizeSweep(doublet_paramSweep, GT=F)
  doublet_BCmvn <- find.pK(doublet_sweep.stats)
  homotypic.prop <- modelHomotypic(GE.data[[i]]@meta.data$seurat_clusters)
  nExp_poi <- round(expectedDoubletPercent[[i]]*length(colnames(GE.data[[i]])))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  print(doublet_BCmvn)
  pK <- as.numeric(as.character(doublet_BCmvn[["pK"]][[which.max(doublet_BCmvn[["BCmetric"]])]]))
  print(pK)
  
  GE.data[[i]]  <- doubletFinder_v3(GE.data[[i]], PCs=1:dims_PCs[[i]], pN=pN, pK=pK, nExp=nExp_poi, reuse.pANN=FALSE, sct=T)
  head(GE.data[[i]]@meta.data)
  classification <- paste("DF.classifications", pN, pK, nExp_poi, sep="_")  # possible error if two have the same classification
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", sample.names[[i]], ".DimPlotClassified.pdf", sep=""),
                   to_plot=DimPlot(GE.data[[i]], group.by=classification))
  print(GE.data[[i]])
  GE.data[[i]] <- SetIdent(GE.data[[i]], value=classification)
  GE.data[[i]] <- subset(GE.data[[i]], idents="Singlet")
  print(GE.data[[i]])
}
saveRDS(GE.data, file=snakemake@output[[1]])