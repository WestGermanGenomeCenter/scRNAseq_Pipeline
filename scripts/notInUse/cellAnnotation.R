library(Seurat)
library(SingleR)

wd <- snakemake@params[[1]]
label.name <- snakemake@params[[2]]
missing.label <- snakemake@params[[3]]
chosen.method <- snakemake@params[[4]]
reference <- readRDS("../", snakemake@input[[1]], sep="")
GE.data <- readRDS("../", snakemake@input[[2]], sep="")

#prepare test reference
reference <- reference[,!is.na(reference[[label.name]]) & reference[[label.name]]!=missing.label]

testset <- as.SingleCellExperiment(GE.data)
results <- SingleR(test=GE.testset, ref=reference, labels=reference[[label.name]], de.method="wilcox") #de.method=chosen.method)
GE.data[["SingleR.labels"]] <- results$labels

saveRDS(GE.data, file=paste("../", snakemake@output[[1]], sep=""))