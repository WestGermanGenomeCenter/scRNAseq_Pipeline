library("Seurat")
#library(future)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

GE <- readRDS(snakemake@input[[1]])
project.name <- snakemake@params[[2]]
#cores <- snakemake@params[[3]]
projectDirPath <- snakemake@params[[1]]

print(Idents(GE))
print(length(GE))

DefaultAssay(GE) <- "RNA"
#plan("multiprocess", workers=cores)
GE <- NormalizeData(GE, assay="RNA")
#plan("multiprocess", workers=cores)
GE <- ScaleData(GE, assay="RNA")
print(GE)
GE.markers <- FindAllMarkers(GE, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.csv(GE.markers, file=paste(projectDirPath, "csv/", project.name,".compareGE_inAll.markers.csv", sep=""))
GE.avgExpression <- AverageExpression(GE,assays="RNA")
write.csv(GE.avgExpression, file=paste(projectDirPath, "csv/", project.name, ".compareGE_inAll.avgExpression.csv", sep=""))

saveRDS(GE, file=snakemake@output[[1]])
