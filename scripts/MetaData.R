library(Seurat)
source("scripts/helperFunctions.R")

raw.data <- Read10X(paste(snakemake@input[[1]], sep=""))
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))
sample.names <- snakemake@params[[2]]
print(sample.names)
project.name <- snakemake@params[[3]]

# create Seurat object and add metadata
if("Gene Expression" %in% names(raw.data)) {
  GE <- CreateSeuratObject(counts=raw.data$`Gene Expression`, project=project.name, min.cells=3, min.features=200, names.field=2, names.delim="-")
  saveRDS(raw.data, file=snakemake@output[[2]]) #file=snakemake@output[[length(snakemake@output)]])
} else {
  GE <- CreateSeuratObject(counts=raw.data, project=project.name, min.cells=3, min.features=200, names.field=2, names.delim="-")
}
head(GE@meta.data)
tail(GE@meta.data)

# divide into 6 original samples
GE.data = SplitObject(GE, split.by="ident")
for(i in 1:length(GE.data)) {
  GE.data[[i]] <- AddMetaData(GE.data[[i]], metadata=sample.names[[i]], col.name="sample")
  saveRDS(GE.data[[i]], file=snakemake@output[[i]])
}
#saveRDS(GE.data, file=snakemake@output[[1]])
