library(Seurat)
library(DoubletFinder)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))
GE.data <- readRDS(snakemake@input[[1]])
sample.name <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
expectedDoubletPercent <- snakemake@params[[3]]
dims_PC <- snakemake@params[[4]]
pN <- 0.25
  

GE.data <- RunUMAP(GE.data, dims=1:dims_PC, verbose=FALSE)
doublet_paramSweep <- paramSweep_v3(GE.data, PCs=1:dims_PC, sct=TRUE)
doublet_sweep.stats <- summarizeSweep(doublet_paramSweep, GT=F)
doublet_BCmvn <- find.pK(doublet_sweep.stats)
homotypic.prop <- modelHomotypic(GE.data@meta.data$seurat_clusters)
nExp_poi <- round(expectedDoubletPercent*length(colnames(GE.data)))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
print(doublet_BCmvn)
pK <- as.numeric(as.character(doublet_BCmvn[["pK"]][[which.max(doublet_BCmvn[["BCmetric"]])]]))
print(pK)
  
GE.data  <- doubletFinder_v3(GE.data, PCs=1:dims_PC, pN=pN, pK=pK, nExp=nExp_poi, reuse.pANN=FALSE, sct=T)
head(GE.data@meta.data)
classification <- paste("DF.classifications", pN, pK, nExp_poi, sep="_")  # possible error if two have the same classification
plot_in_terminal(plotname=paste(projectDirPath, "plots/", sample.name, ".DimPlotClassified.pdf", sep=""),
                 to_plot=DimPlot(GE.data, group.by=classification))
print(GE.data)
GE.data <- SetIdent(GE.data, value=classification)
GE.data <- subset(GE.data, idents="Singlet")
print(GE.data)
saveRDS(GE.data, file=snakemake@output[[1]])