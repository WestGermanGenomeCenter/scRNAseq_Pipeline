plot_in_terminal <- function(plotname, to_plot) {
  pdf(plotname)
  plot(to_plot)
  dev.off()
  cat(paste("\n", "To check the Plot, see '", plotname, "'\n\n", sep=""))
}