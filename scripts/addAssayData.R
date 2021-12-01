library("Seurat")
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

GE.215.integrated <- readRDS(snakemake@input[[1]])
raw.data <- readRDS(snakemake@input[[2]])
project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
assayname <- snakemake@params[[3]]

#print(colnames(GE.215.integrated))

GE.215.integrated[[assayname]] <- CreateAssayObject(counts=raw.data$`Antibody Capture`[, colnames(GE.215.integrated)])
GE.215.integrated <- NormalizeData(GE.215.integrated, assay=assayname, normalization.method="CLR")
GE.215.integrated <- ScaleData(GE.215.integrated, assay=assayname)

saveRDS(GE.215.integrated, file=snakemake@output[[1]])