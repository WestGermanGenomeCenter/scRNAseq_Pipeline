library("Seurat")
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))

GE.integrated <- readRDS(snakemake@input[[1]])
raw.data <- readRDS(snakemake@input[[2]])
project.name <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
assay.name <- snakemake@params[[3]]

GE.integrated[[assay.name]] <- CreateAssayObject(counts=raw.data$`Antibody Capture`[, colnames(GE.integrated)])
GE.integrated <- NormalizeData(GE.integrated, assay=assay.name, normalization.method="CLR")
GE.integrated <- ScaleData(GE.integrated, assay=assay.name)

saveRDS(GE.integrated, file=snakemake@output[[1]])