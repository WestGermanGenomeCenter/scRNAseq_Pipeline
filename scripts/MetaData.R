library(Seurat)

raw.data <- Read10X(paste(snakemake@input[[1]], sep=""))
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))
names <- snakemake@params[[2]]
project <- snakemake@params[[3]]

GE215 <- 0
# create Seurat object and add metadata
if("Gene Expression" %in% names(raw.data)) {
  GE215 <- CreateSeuratObject(counts=raw.data$`Gene Expression`, project=project, min.cells=3, min.features=200, names.field=2, names.delim="-")
  saveRDS(raw.data, file=snakemake@output[[2]])
} else {
  GE215 <- CreateSeuratObject(counts=raw.data, project=project, min.cells=3, min.features=200, names.field=2, names.delim="-")
}
head(GE215@meta.data)
tail(GE215@meta.data)

# divide into 6 original samples
GE_Data = SplitObject(GE215, split.by="ident")
for(i in 1:length(GE_Data)) {
  GE_Data[[i]] <- AddMetaData(GE_Data[[i]], metadata = names[[i]], col.name = "sample")
}
saveRDS(GE_Data, file=snakemake@output[[1]])
#saveRDS(GE215, file=paste("../", snakemake@output[[2]], sep=""))
