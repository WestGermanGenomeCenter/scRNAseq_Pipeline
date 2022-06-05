import networkx as nx
import pygraphviz as pyviz
import matplotlib.pyplot as plt
import plotters.installRunPlot as iplot
import plotters.pipelineRunPlot as pPlot

instRun = nx.nx_agraph.to_agraph(iplot.installRun)
instRun.layout("dot")
instRun.draw("figures/installRun.svg")
pipeRun = nx.nx_agraph.to_agraph(pPlot.pipelineRun)
pipeRun.layout("dot")
pipeRun.draw("figures/pipelineRun.svg")

