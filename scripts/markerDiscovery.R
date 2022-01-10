library("Seurat")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

GE.215.integrated <- readRDS(snakemake@input[[1]])
project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]

print(Idents(GE.215.integrated))
print(length(GE.215.integrated))

DefaultAssay(GE.215.integrated) <- "RNA"
GE.215.integrated <- NormalizeData(GE.215.integrated, assay="RNA")
GE.215.integrated <- ScaleData(GE.215.integrated, assay="RNA")
print(GE.215.integrated)
GE.215.integrated.markers <- FindAllMarkers(GE.215.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(GE.215.integrated.markers, file=paste(projectDirPath, "csv/", project,".compareGE_inAll.markers.csv", sep=""))
GE.215.integrated.avgExpression <- AverageExpression(GE.215.integrated,assays = "RNA")
write.csv(GE.215.integrated.avgExpression, file=paste(projectDirPath, "csv/", project, ".compareGE_inAll.avgExpression.csv", sep=""))

saveRDS(GE.215.integrated, file=snakemake@output[[1]])
