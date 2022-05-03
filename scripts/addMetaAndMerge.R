library(Seurat)
library(gtools)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))

GE.data <- snakemake@input
#GE <- readRDS(paste("../", snakemake@input[[2]], sep=""))
condition <- snakemake@params[[2]]
condition.names <- snakemake@params[[3]]
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

num.conditions <- length(condition.names)
if(num.conditions > 1) {
  for(i in 2:num.conditions) {
    perms <- permutations(num.conditions, i, condition.names)
    print(perms)
    for(j in 1:nrow(perms)) {
      split.list <- SplitObject(GE, split.by=perms[[j]][[1]])
      for(k in 2:ncol(perms)) {
        for(l in 1:length(split.list)) {
          tmp.list[[l]] <- (split.list[l], split.by=perms[[j]][[k]])
        }
        split.list <- unlist(tmp.list, recursive=TRUE)
      }
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
