plot_in_terminal <- function(plotname, to_plot) {
  pdf(plotname)
  plot(to_plot)
  dev.off()
  cat(paste("\n", "To check the Plot, see '", plotname, "'\n\n", sep=""))
}

#https://stackoverflow.com/questions/64101921/snakemake-how-to-log-python-scripts-executed-by-script-directive
#redirect, so snakemake log and script directive can be used
log <- file(snakemake@log[[1]], open="wt")
sink(log)