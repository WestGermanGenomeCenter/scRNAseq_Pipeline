import networkx as nx
import pygraphviz as pyviz
import matplotlib.pyplot as plt
import plotters.parallelBarplot as bPlot
import plotters.installRunPlot as iPlot
import plotters.pipelineRunPlot as pPlot

instRun = nx.nx_agraph.to_agraph(iPlot.installRun)
instRun.layout("dot")
instRun.draw("figures/installRun.pdf")
pipeRun = nx.nx_agraph.to_agraph(pPlot.pipelineRun)
pipeRun.layout("fdp")
pipeRun.draw("figures/pipelineRun.pdf")

bPlot.data.pivot("rule", "sample_numCores", "time_in_s").plot(kind="bar")
plt.xticks(rotation=0)
plt.ylabel("walltime in seconds")
plt.savefig("figures/paraFutureBars.png")
