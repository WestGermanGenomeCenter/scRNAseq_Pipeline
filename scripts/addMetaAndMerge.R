library(Seurat)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))

GE.data <- readRDS(snakemake@input[[1]])
#GE215 <- readRDS(paste("../", snakemake@input[[2]], sep=""))
condition <- snakemake@params[[2]]
condition.name <- snakemake@params[[3]]

for(i in 1:length(GE.data)) {
  GE.data[[i]] <- AddMetaData(GE.data[[i]], metadata=condition[[i]], col.name=condition.name)
}

print(head(GE.data[[1]]@meta.data))
#print(GE215)
if(length(GE.data) > 1) {
  GE215 <- merge(GE.data[[1]], y = GE.data[2:length(GE.data)])
} else {
  GE215 <- GE.data[[1]]
}
DefaultAssay(GE215) <- "RNA"
print(GE215)
print(head(GE215@meta.data))

saveRDS(GE215, file=snakemake@output[[1]])
