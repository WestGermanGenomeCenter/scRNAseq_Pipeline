library(Seurat)
library(DoubletFinder)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))
GE_Data <- readRDS(snakemake@input[[1]])
names <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
roundingValues <- snakemake@params[[3]]
dims_PCs <- snakemake@params[[4]]
pN <- 0.25
  
for(i in 1:length(GE_Data)) {
  GE_Data[[i]] <- RunUMAP(GE_Data[[i]], dims=1:dims_PCs[[i]], verbose=FALSE)
#  plot_in_terminal(plotname=paste(projectDirPath, "plots/", names[[i]], ".DimPlot.pdf", sep=""),
#                   to_plot=DimPlot(GE_Data[[i]]))
  doublet_paramSweep <- paramSweep_v3(GE_Data[[i]], PCs=1:dims_PCs[[i]], sct=TRUE)
  doublet_sweep.stats <- summarizeSweep(doublet_paramSweep, GT=F)
  doublet_BCmvn <- find.pK(doublet_sweep.stats)
  homotypic.prop <- modelHomotypic(GE_Data[[i]]@meta.data$seurat_clusters)
  nExp_poi <- round(roundingValues[[i]]*length(colnames(GE_Data[[i]])))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  print(doublet_BCmvn)
  pK <- as.numeric(as.character(doublet_BCmvn[["pK"]][[which.max(doublet_BCmvn[["BCmetric"]])]]))
  print(pK)
  
  GE_Data[[i]]  <- doubletFinder_v3(GE_Data[[i]], PCs=1:dims_PCs[[i]], pN=pN, pK=pK, nExp=nExp_poi, reuse.pANN=FALSE, sct=T)
  head(GE_Data[[i]]@meta.data)
  classification <- paste("DF.classifications", pN, pK, nExp_poi, sep="_")  # possible error if two have the same classification
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", names[[i]], ".DimPlotClassified.pdf", sep=""),
                   to_plot=DimPlot(GE_Data[[i]], group.by=classification))
  print(GE_Data[[i]])
  GE_Data[[i]] <- SetIdent(GE_Data[[i]], value=classification)
  GE_Data[[i]] <- subset(GE_Data[[i]], idents="Singlet")
#  plot_in_terminal(plotname=paste(projectDirPath, "plots/", names[[i]], ".DimPlotClassifiedSinglet.pdf", sep=""),
#                   to_plot=DimPlot(GE_Data[[i]], group.by=classification))
  print(GE_Data[[i]])
}
saveRDS(GE_Data, file=snakemake@output[[1]])