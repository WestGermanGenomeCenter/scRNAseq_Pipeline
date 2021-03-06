library("Seurat")
library(ggplot2)
library(seurathelpeR)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))


GE <- readRDS(snakemake@input[[1]])
condition <- snakemake@params[[2]]
projectDirPath <- snakemake@params[[1]]
res_ident <- paste("integrated_snn_res", snakemake@params[[3]], sep=".")
project <- snakemake@params[[5]]
freq <- "freq"
sdplus <- "sdplus"
sdminus <- "sdminus"

if(snakemake@params[[4]]) {
  ccount <- count_cells(seurat_object=GE, group_by_var=condition, subgroup_var=res_ident)
  ccount[["freq"]] <- ccount[["perc"]]/100
  
  #Errorbars
  GE <- SplitObject(GE, split.by=condition)
  df <- data.frame(ecount=numeric(0), condition.name=character(0), cluster.num=numeric(0))
  for(i in 1:length(GE)) { #split by condition
    e <- count_cells(seurat_object=GE[[i]], group_by_var="orig.ident", subgroup_var=res_ident)
    condition.name <- GE[[i]]@meta.data[[condition]][[1]]
    for(j in 1:length(unique(e[[res_ident]]))) {
      filteredE <- e[e[, res_ident] == j-1, ]$perc
      if(length(filteredE) < length(unique(e$orig.ident))) { #if some samples are not contained in clusters add zeroes
        filteredE <- c(filteredE, numeric(length(unique(e$orig.ident))-length(filteredE)) )
      }
      ecount <- sd(filteredE/100)
      df <- rbind(df, data.frame(res_ident=j-1, condition.name, ecount))
    }
  }
  colnames(df) <- c(res_ident, condition, "ecount")
  ccount <- merge(df, ccount, by=c(res_ident, condition))
  ccount[[sdplus]] <- ccount[[freq]] + ccount[["ecount"]]
  ccount[[sdminus]] <- ccount[[freq]] - ccount[["ecount"]]
  
  #conditions together
  plot <- ggplot(ccount, aes_string(x=res_ident, y=freq, fill=condition)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(position=position_dodge(), aes_string(x=res_ident, ymin=sdminus, ymax=sdplus))
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".", condition, ".barplot.pdf", sep=""),
                   to_plot=plot)
  
  #individual conditions
  conditions <- unique(ccount[[condition]])
  for(i in 1:length(conditions)) {
    plot <- ggplot(ccount[ccount[, condition] == conditions[[i]],], aes_string(x=res_ident, y=freq)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(position=position_dodge(), aes_string(x=res_ident, ymin=sdminus, ymax=sdplus))
    plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".", conditions[[i]], ".barplot.pdf", sep=""),
                     to_plot=plot)
  }
} else {
  res_ident <- "seurat_clusters"
  ccount <- count_cells(seurat_object=GE, group_by_var=res_ident)
  ccount[["freq"]] <- ccount[["perc"]]/100
  plot <- ggplot(ccount, aes_string(x=res_ident, y=freq)) + geom_bar(position=position_dodge(), stat="identity")
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".", condition, ".barplot.pdf", sep=""),
                   to_plot=plot)
}
write.csv(ccount, file=paste(projectDirPath, "csv/", project, ".cellCountingTable.csv", sep=""))