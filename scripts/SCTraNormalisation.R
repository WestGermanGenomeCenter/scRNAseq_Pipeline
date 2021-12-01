library(Seurat)
library(sctransform)

setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))
GE <- readRDS(snakemake@input[[1]])
GE.215.list <- SplitObject(GE, split.by = "sample")

for (i in 1:length(GE.215.list)) {
  GE.215.list[[i]] <- SCTransform(GE.215.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}

saveRDS(GE.215.list, file=snakemake@output[[1]])