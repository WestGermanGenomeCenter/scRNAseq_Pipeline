library(Seurat)
library(sctransform)
source("scripts/helperFunctions.R")

setwd(paste(snakemake@params[[1]], "workDirectory", sep=""))
GE <- readRDS(snakemake@input[[1]])
GE.data <- SplitObject(GE, split.by="sample")

for (i in 1:length(GE.data)) {
  GE.data[[i]] <- SCTransform(GE.data[[i]], vars.to.regress="percent.mt", verbose=FALSE)
}

saveRDS(GE.data, file=snakemake@output[[1]])