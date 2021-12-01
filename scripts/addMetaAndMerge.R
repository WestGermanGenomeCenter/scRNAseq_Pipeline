library(Seurat)
setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))

GE_Data <- readRDS(snakemake@input[[1]])
#GE215 <- readRDS(paste("../", snakemake@input[[2]], sep=""))
meta <- snakemake@params[[2]]
colname <- snakemake@params[[3]]

for(i in 1:length(GE_Data)) {
  GE_Data[[i]] <- AddMetaData(GE_Data[[i]], metadata=meta[[i]], col.name=colname)
}

print(head(GE_Data[[1]]@meta.data))
#print(GE215)
if(length(GE_Data) > 1) {
  GE215 <- merge(GE_Data[[1]], y = GE_Data[2:length(GE_Data)])
} else {
  GE215 <- GE_Data[[1]]
}
DefaultAssay(GE215) <- "RNA"
print(GE215)
print(head(GE215@meta.data))

saveRDS(GE215, file=snakemake@output[[1]])
