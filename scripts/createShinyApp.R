library(Seurat)
library(ShinyCell)
source("scripts/helperFunctions.R")

project <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]

if(snakemake@params[[3]]) {
  gexAssay = "integrated"
} else {
  gexAssay <- "RNA"
}

setwd(snakemake@params[[1]])
seu <- readRDS(snakemake@input[[1]])
pathToShinyApp <- paste(projectDirPath, "shinyApp/", sep="")
scConf <- createConfig(seu)

metaToDel <- c()
j <- 1
for(i in 1:length(scConf$ID)) {
  if(grepl("pANN", scConf$ID[[i]], fixed=TRUE)) {
    metaToDel[[j]] <- scConf$ID[[i]]
    j <- j + 1
  }
  if(grepl("DF.classifications", scConf$ID[[i]], fixed=TRUE)) {
    metaToDel[[j]] <- scConf$ID[[i]]
    j <- j + 1
  }
}

scConf <- delMeta(scConf, metaToDel)
makeShinyApp(seu, scConf, gene.mapping=T, gex.assay=gexAssay,
             shiny.title=paste("ShinyCell", project),
             shiny.dir=pathToShinyApp)