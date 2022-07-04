library("Seurat")
source("scripts/helperFunctions.R")
GE <- Read10X(snakemake@input[[1]])

setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))
project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
assayname <- snakemake@params[[3]]
names <- snakemake@params[[4]]

HTO <- GE$`Antibody Capture`
GE <- GE$`Gene Expression`

joint.bcs <- intersect(colnames(GE), colnames(HTO))
GE <- GE[, joint.bcs]
HTO <- as.matrix(HTO[, joint.bcs])
#rows <- rownames(HTO)
#rows

GE <- CreateSeuratObject(counts=GE)
GE <- SCTransform(GE, verbose=FALSE)
GE[[assayname]] <- CreateAssayObject(counts=HTO)
GE <- NormalizeData(GE, assay=assayname, normalization.method="CLR")

GE <- HTODemux(GE, assay=assayname, positive.quantile=0.99)

table(GE$HTO_classification.global)
Idents(GE) <- "HTO_maxID"
rows <- levels(Idents(GE))
#FeatureScatter(GE, feature1=paste("hto", rows[[1]], sep="_"), feature2=paste("hto", rows[[2]], sep="_"))
plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".HTOscatter.pdf", sep=""),
                 to_plot=FeatureScatter(GE, feature1=paste("hto", rows[[1]], sep="_"), feature2=paste("hto", rows[[2]], sep="_")))

Idents(GE) <- "HTO_classification.global"
GE <- subset(GE, idents="Doublet", invert=TRUE)

Idents(GE) <- "HTO_maxID"
GE <- SplitObject(GE)

for(i in 1:length(GE)) {
  GE[[i]] <- AddMetaData(GE[[i]], metadata = names[[i]], col.name = "sample")
  saveRDS(GE[[i]], file=paste(snakemake@output[[i]]))
}

#saveRDS(GE, snakemake@output[[1]])
