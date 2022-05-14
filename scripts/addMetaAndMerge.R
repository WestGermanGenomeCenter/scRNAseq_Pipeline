library(Seurat)
library(gtools)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))

GE.data <- snakemake@input
#GE <- readRDS(paste("../", snakemake@input[[2]], sep=""))
condition <- snakemake@params[[2]]
condition.names <- snakemake@params[[3]]
condition.combis <- snakemake@params[[4]]
print(condition)
print(condition.names)
sampleNum <- 0

for(i in 1:length(GE.data)) {
  GE.data[[i]] <- readRDS(GE.data[[i]])
  for(j in 1:length(condition.names)) {
    print(condition[[sampleNum + j]])
    print(condition.names[[j]])
    GE.data[[i]] <- AddMetaData(GE.data[[i]], metadata=condition[[sampleNum + j]], col.name=condition.names[[j]])
  }
  sampleNum <- sampleNum + length(condition.names)
}

print("\nNow combined meta")
if(length(condition.names) > 1) {
  print(condition.combis)
  for(j in 1:length(GE.data)) {
    for(k in 1:length(condition.combis)) {
      new.meta <- condition.combis[[k]]
      print(new.meta)
      meta <- unlist(strsplit(new.meta, "_"))
      print(meta)
      for(l in 1:length(meta)) {
        meta[[l]] <- GE.data[[j]]@meta.data[[meta[[l]]]][[1]]
	print(meta[[l]])
      }
      meta <- paste(meta, collapse="_")
      print(meta)
      GE.data[[j]] <- AddMetaData(GE.data[[j]], metadata=meta, col.name=new.meta)
    }
  }
}



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
