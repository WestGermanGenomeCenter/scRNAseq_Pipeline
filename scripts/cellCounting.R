library("Seurat")
library(ggplot2)
source("scripts/helperFunctions.R")
setwd(paste(snakemake@params[[1]], "workDirectory/", sep=""))


GE.integrated <- readRDS(snakemake@input[[1]])
condition <- snakemake@params[[2]]
Idents(GE.integrated) <- condition
projectDirPath <- snakemake@params[[1]]
res_ident <- paste("integrated_snn_res", snakemake@params[[3]], sep=".")
project <- snakemake@params[[5]]
ids <- levels(Idents(GE.integrated))
Idents(GE.integrated) <- res_ident
num_clusters <- length(levels(Idents(GE.integrated)))
comb_bars <- data.frame()

if(snakemake@params[[4]]) {
  GE.integrated <- SplitObject(GE.integrated, split.by=condition)
  for(i in 1:length(ids)) {
    #bars
    bars <- data.frame(table(Idents(GE.integrated[[i]])))
    bars$Freq <- bars$Freq/sum(bars$Freq)
    names(bars)[names(bars) == "Var1"] <- "cluster"
    while(length(bars$Freq) < num_clusters) {
      tmp <- length(bars)-1
      new_last_row <- data.frame(Var1=as.factor(tmp), Freq=0)
      names(new_last_row) <- c("cluster", "Freq")
      bars <- rbind(bars, new_last_row)
    }
    #error bars
    GE.tables <- SplitObject(GE.integrated[[i]], split.by="orig.ident")
    for(j in 1:length(GE.tables)) { 
      #Split into samples belonging to condtition and turn them relative Freq
      GE.tables[[j]] <- data.frame(table(Idents(GE.tables[[j]])))
      GE.tables[[j]]$Freq <- GE.tables[[j]]$Freq/sum(GE.tables[[j]]$Freq)
      len_table <- length(GE.tables[[j]]$Freq)
      #if the current sample does not have any cells in the last cluster
      if(len_table < num_clusters & GE.tables[[j]]$Var1[[len_table]] == len_table-1) {
        tmp <- num_clusters-1
        new_last_row <- data.frame(Var1=as.factor(tmp), Freq=0)
        names(new_last_row) <- c("Var1", "Freq")
        GE.tables[[j]] <- rbind(GE.tables[[j]], new_last_row)
      }
    }
    s_deviation <- c(1:num_clusters)
    for(j in 1:num_clusters) {
      curr_clusters <- c(1:length(GE.tables))
      #fill in for replicates with no cells in one cluster
      for(k in 1:length(GE.tables)) {
        if(GE.tables[[k]]$Var1[[j]] == j-1) {
          curr_clusters[[k]] <- GE.tables[[k]]$Freq[[j]]
        } else {
          GE.tables <- rbind(GE.tables[1:j-1, ], GE.tables[1:1, ]*0, GE.tables[j:num_clusters])
          curr_clusters[[k]] <- GE.tables[[k]]$Freq[[j]]
        }
        s_deviation[[j]] <- sd(curr_clusters)
      }
    }
    #plot
    plot <- ggplot(bars) + geom_bar(aes(x=cluster, y=Freq), stat="identity") + geom_errorbar(aes(x=cluster, ymin=Freq-s_deviation, ymax=Freq+s_deviation))
    plot_in_terminal(plotname=paste(projectDirPath, "plots/", ids[[i]], ".barplot.pdf", sep=""),
                     to_plot=plot)
    print(plot)
    bars[["sd"]] <- s_deviation
    bars[["condition"]] <- rep(ids[[i]], num_clusters)
    comb_bars <- rbind(comb_bars, bars)
  }
  plot <- ggplot(comb_bars, aes(x=cluster, y=Freq, fill=condition)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(position=position_dodge(), aes(x=cluster, ymin=Freq-sd, ymax=Freq+sd))
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".", condition, ".barplot.pdf", sep=""),
                   to_plot=plot)
} else {
  comb_bars <- data.frame(table(Idents(GE.integrated)))
  comb_bars$Freq <- comb_bars$Freq/sum(comb_bars$Freq)
  names(comb_bars)[names(comb_bars) == "Var1"] <- "cluster"
  plot <- ggplot(comb_bars, aes(x=cluster, y=Freq, fill=condition)) + geom_bar(position=position_dodge(), stat="identity")
  plot_in_terminal(plotname=paste(projectDirPath, "plots/", project, ".", condition, ".barplot.pdf", sep=""),
                   to_plot=plot)
}

