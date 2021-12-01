# Preprocessing (addition of Metadata & QC)
setwd("data\\Seurat")
library(Seurat)
library(DoubletFinder)
raw.data <- Read10X("..\\215_cellranger_aggr\\outs\\filtered_feature_bc_matrix")
# create Seurat object and add metadata
GE215 <- CreateSeuratObject(counts = raw.data, project = "215", min.cells = 3, min.features = 200,names.field = 2,names.delim = "-")
head(GE215@meta.data)
tail(GE215@meta.data)
#divide into 6 original samples
                        names <- list("157-1", "215-1", "215-2", "215-3", "215-4","215-5")
GE_Data = SplitObject(GE215, split.by="ident")
for(i in 1:length(GE_Data)) {
  GE_Data[[i]] <- AddMetaData(GE_Data[[i]], metadata = names[[i]], col.name = "sample")
}

#Mitochondrial Analysis
for(i in 1:length(GE_Data)) {
  GE_Data[[i]][["percent.mt"]] <- PercentageFeatureSet(GE_Data[[i]], pattern = "^mt-")
  print(VlnPlot(GE_Data[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  print(mean(GE_Data[[i]]@meta.data$percent.mt))
  print(median(GE_Data[[i]]@meta.data$percent.mt))
  cutoff <- as.numeric(readline("Input the cutoff for percent.mt based on mean and median: "))
  #cat("Input the cutoff for percent.mt based on mean and median: ")
  #cutoff <- as.numeric(readLines("stdin", 1))
  print(names[[i]])
  print(GE_Data[[i]])
  #GE_Data[[i]] <- subset(GE_Data[[i]], subset = percent.mt < 10)
  GE_Data[[i]] <- subset(GE_Data[[i]], subset = percent.mt < cutoff)
  print(GE_Data[[i]])
  print(VlnPlot(GE_Data[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  readline("Time to compare Plots, press ENTER to continue")
}
save.image("215_Seurat_Analysis.RData")
savehistory("215_Seurat_Analysis.Rhistory")

#Doublet Removal
                        pN <- 0.25
                        pKs <- list(0.005, 0.08, 0.09, 0.3, 0.09, 0.005)
                        rounding <- list(0.02, 0.01, 0.016, 0.016, 0.016, 0.016)
#! nExp_poi not correct for nExp? Or because of 1:17 instead of 1:30
for(i in 1:length(GE_Data)) {
  GE_Data[[i]] <- SCTransform(GE_Data[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
  GE_Data[[i]] <- RunPCA(GE_Data[[i]], verbose = FALSE)
  print(ElbowPlot(GE_Data[[i]]))
  GE_Data[[i]] <- RunUMAP(GE_Data[[i]], dims = 1:17, verbose = FALSE)
  print(DimPlot(GE_Data[[i]]))
  doublet_paramSweep <- paramSweep_v3(GE_Data[[i]], PCs = 1:30, sct = TRUE) #here 1:17 should actually be used
  doublet_sweep.stats <- summarizeSweep(doublet_paramSweep, GT=F)
  doublet_BCmvn <- find.pK(doublet_sweep.stats)
  homotypic.prop <- modelHomotypic(GE_Data[[i]]@meta.data$seurat_clusters)
  nExp_poi <- round(rounding[[i]]*length(colnames(GE_Data[[i]])))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  View(doublet_BCmvn, title=paste(names[[i]], "Doublet_BCmvn", sep="_"))
  #pK <- as.numeric(readline("Enter pK: "))
  #GE_Data[[i]]  <- doubletFinder_v3(GE_Data[[i]], PCs = 1:30, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  GE_Data[[i]]  <- doubletFinder_v3(GE_Data[[i]], PCs = 1:30, pN = pN, pK = pKs[[i]], nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  head(GE_Data[[i]]@meta.data)
  #classification <- paste("DF.classifications", pN, pK, nExp_poi, sep="_")
  classification <- paste("DF.classifications", pN, pKs[[i]], nExp_poi, sep="_")
  print(DimPlot(GE_Data[[i]],group.by = classification))
  GE_Data[[i]]
  GE_Data[[i]] <- SetIdent(GE_Data[[i]],value = classification)
  GE_Data[[i]] <- subset(GE_Data[[i]],idents = "Singlet")
  print(DimPlot(GE_Data[[i]],group.by = classification))
  GE_Data[[i]]
}

#Add timepoints
                        timepoints <- list("Tag3", "Tag3", "Tag3", "Tag7", "Tag7", "Tag7")
for(i in 1:length(GE_Data)) {
  GE_Data[[i]] <- AddMetaData(GE_Data[[i]], metadata=timepoints[[i]], col.name="timepoint")
}

#merge back together, make it more generalized
print(head(GE_Data[[1]]@meta.data))
GE215
GE215 <- merge(GE_Data[[1]], y = c(GE_Data[[2]], GE_Data[[3]], GE_Data[[4]], GE_Data[[5]], GE_Data[[6]]))
GE215
DefaultAssay(GE215) <- "RNA"
GE215
head(GE215@meta.data)
save.image("215_Seurat_preprocessing.RData")
save(GE215, file = "GE.215.preprocessed.Robj")

###############################################Below untested

#LOAD OBJECT IN NEW R SESSION
#Normalisation via scTransfrom
load("GE.215.preprocessed.Robj")
library(sctransform)
GE.215.list <- SplitObject(GE215, split.by = "sample")
for (i in 1:length(GE.215.list)) {
  GE.215.list[[i]] <- SCTransform(GE.215.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}

#Cannot test below due to not having enough RAM
# Integration and dimensionality reduction
GE.215.features <- SelectIntegrationFeatures(GE.215.list, nfeatures = 3000)
options(future.globals.maxSize= 4294967296)
GE.215.list <- PrepSCTIntegration(object.list = GE.215.list, anchor.features = GE.215.features, verbose = FALSE)
GE.215.anchors <- FindIntegrationAnchors(object.list = GE.215.list, normalization.method = "SCT", anchor.features = GE.215.features, verbose = FALSE)
GE.215.integrated <- IntegrateData(anchorset = GE.215.anchors, normalization.method = "SCT", verbose = FALSE) #memory error
GE.215.integrated <- RunPCA(GE.215.integrated, verbose = FALSE)
GE.215.integrated <- RunUMAP(GE.215.integrated, dims = 1:30)
DimPlot(GE.215.integrated,group.by = "sample")
ElbowPlot(GE.215.integrated,ndims = 30)
GE.215.integrated <- FindNeighbors(GE.215.integrated, dims = 1:25)
continue <- TRUE
while(continue) {
  resolution <- as.numeric(readline("Input desired resolution: "))
  GE.215.integrated <- FindClusters(GE.215.integrated, resolution=resolution)
  print(DimPlot(GE.215.integrated))
  happy <- readline("Are you happy with the resolution? [y/n]")
  if(happy == "y") {
    continue <- FALSE
  }
}
DimPlot(GE.215.integrated, label=T)
DimPlot(GE.215.integrated,split.by = "timepoint")

#Marker Discovery
DefaultAssay(GE.215.integrated) <- "RNA"
GE.215.integrated <- NormalizeData(GE.215.integrated,assay = "RNA")
GE.215.integrated <- ScaleData(GE.215.integrated,assay = "RNA")
GE.215.integrated
GE.215.integrated.markers <- FindAllMarkers(GE.215.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(GE.215.integrated.markers,file = "./GE.215.integrated.markers.csv")
GE.215.integrated.avgExpression <- AverageExpression(GE.215.integrated,assays = "RNA")
write.csv(GE.215.integrated.avgExpression,file = "GE.215.integrated.avgExpression.csv")

#Cell counting
Idents(GE.215.integrated) <- "orig.ident"
for(i in 1:6) {
  GE.integrated <- subset(GE.215.integrated, idents=i)
  Idents(GE.integrated) <- "integrated_snn_res.0.8"
  GE.integrated.cellCount = table(Idents(GE.integrated))
  write.csv(GE.integrated.cellCount, file=paste(names[[i]], "_cellCount.csv"))
}

                    timepoints2 <- list("Tag3", "Tag7")
Idents(GE.215.integrated) <- "timepoint"
for(i in 1:length(timepoints2)) {
  GE.integrated.timepoint <- subset(GE.215.integrated, idents=timepoints2[[i]])
  print(GE.integrated.timepoint)
  Idents(GE.integrated.timepoint) <- "integrated_snn_res.0.8"
  GE.integrated.timepoint.cellCount <- table(Idents(GE.integrated.timepoint))
  write.csv(GE.integrated.timepoint.cellCount, file = paste("GE.integrated", timepoints2[[i]], ".csv"))
}
Idents(GE.215.integrated) <- "integrated_snn_res.0.8"
save.image("215_Seurat_Analysis.RData")