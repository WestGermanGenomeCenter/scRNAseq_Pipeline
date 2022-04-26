library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))

GE.data <- snakemake@input
#GE <- readRDS(paste("../", snakemake@input[[2]], sep=""))
condition <- snakemake@params[[2]]
condition.name <- snakemake@params[[3]]
print(condition)
print(condition.name)
sampleNum <- 0

for(i in 1:length(GE.data)) {
  GE.data[[i]] <- readRDS(GE.data[[i]])
  for(j in 1:length(condition.name)) {
    print(condition[[sampleNum + j]])
    print(condition.name[[j]])
    GE.data[[i]] <- AddMetaData(GE.data[[i]], metadata=condition[[sampleNum + j]], col.name=condition.name[[j]])
  }
  sampleNum <- sampleNum + length(condition.name)
}
print("LOOP DONE!")
print(head(GE.data[[1]]@meta.data))
#print(GE)
if(length(GE.data) > 1) {
  GE <- merge(GE.data[[1]], y = GE.data[2:length(GE.data)])
} else {
  GE <- GE.data[[1]]
}
DefaultAssay(GE) <- "RNA"
print(GE)
print(head(GE@meta.data))

saveRDS(GE, file=snakemake@output[[1]])
