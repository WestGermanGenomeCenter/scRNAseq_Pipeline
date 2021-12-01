library(devtools)
library(withr)
absPath <- paste(snakemake@params[[1]], "/workDirectory", sep="")
doubletPath <- paste("../", snakemake@params[[2]], "/lib/R/library", sep="")
shinyPath <- paste("../", snakemake@params[[3]], "/lib/R/library", sep="")
setwd(absPath)

with_libpaths(new=doubletPath, install_local("../dependencies/DoubletFinder-master.zip"))
with_libpaths(new=shinyPath, install_local("../dependencies/ShinyCell-master.zip"))
