library(devtools)
library(withr)
source("scripts/helperFunctions.R")
absPath <- paste(snakemake@params[[1]], "/workDirectory", sep="")
doubletPath <- paste("../", snakemake@params[[2]], "/lib/R/library", sep="")
shinyPath <- paste("../", snakemake@params[[3]], "/lib/R/library", sep="")
countPath <- paste("../", snakemake@params[[4]], "/lib/R/library", sep="")
hhu_hpc <- snakemake@params[[5]]
setwd(absPath)

with_libpaths(new=doubletPath, install_local("../dependencies/DoubletFinder-master.zip"))
with_libpaths(new=shinyPath, install_local("../dependencies/ShinyCell-master.zip"))

#if not using hhu hpc, these packages should be downloaded via conda from snakemake
#the hhu hpc currently has no internet and issues with some of these packages, rest are needed due to compatibility
if(hhu_hpc) {
    with_libpaths(new=countPath, install_local("../dependencies/GenomeInfoDbData_1.2.7.tar.gz"))
    with_libpaths(new=countPath, install_local("../dependencies/GenomeInfoDb_1.30.0.tar.gz"))
    with_libpaths(new=countPath, install_local("../dependencies/GenomicRanges_1.46.1.tar.gz"))
    with_libpaths(new=countPath, install_local("../dependencies/SummarizedExperiment_1.24.0.tar.gz"))
    with_libpaths(new=countPath, install_local("../dependencies/SingleCellExperiment_1.16.0.tar.gz"))
}
with_libpaths(new=countPath, install_local("../dependencies/seurathelpeR-master.zip"))