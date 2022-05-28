import networkx as nx
import pygraphviz as pyg
import matplotlib.pyplot as plt

#rules = [
#    ("metaData", {"pos": (4,12), "size": 1000}),
#    ("demultiplexing", {"pos": (6,12), "size": 1000}),
#    ("mt_p1", {"pos": (5,11), "size": 1000}),
#    ("mt_p2", {"pos": (5,10), "size": 600}),
#    ("doubletRemovalElbowPlot", {"pos": (5,9), "size": 600}),
#    ("doubletRemoval", {"pos": (5,8), "size": 600}),
#    ("addTPsMerge", {"pos": (5,7), "size": 600}),
#    ("SCTransformNormalization", {"pos": (5,6), "size": 600}),
#    ("IntegrationDimReduction", {"pos": (5,5), "size": 600}),
#    ("RunUMAP", {"pos": (5,4), "size": 1500}),
#    ("testDiffClusterResolutions", {"pos": (3,3), "size": 1000}),
#    ("useChosenClusterResolutions", {"pos": (5, 3), "size": 1000}),
#    ("multimodalAnalysis", {"pos": (8, 2), "size": 1000}),
#    ("markerDiscovery", {"pos": (4, 2), "size": 1000}),
#    ("cellCounting", {"pos": (2,1), "size": 1000}),
#    ("DGE", {"pos": (6,1)}),
#    ("createShinyApp", {"pos": (4,1)}),
#    ("multimodalFeaturePlotting", {"pos": (8,1), "size": 1000})
#]

#installation = nx.DiGraph()
#installation.add_nodes_from([
#    ("createDoubletEnv", {"pos": (1,2)}),
#    ("createCountEnv", {"pos": (2,2)}),
#    ("createShinyEnv", {"pos": (3,2)}),
#    ("installMissingPackages", {"pos": (2,1)})])
#installation.add_edges_from([
#    ("createDoubletEnv", "installMissingPackages"), ("createCountEnv", "installMissingPackages"), ("createShinyEnv", "installMissingPackages")
#])
#pos = nx.get_node_attributes(installation, "pos")
#nx.draw(installation, pos, node_shape="8", with_labels=True)
#plt.show()

pipelineRun = nx.MultiDiGraph()
#pipelineRun.add_nodes_from(rules)
pipelineRun.add_edges_from([
    #ns
    ("metaData", "mt_p1", {"color": "black"}),
    ("mt_p1", "mt_p2", {"color": "black"}),
    ("mt_p2", "doubletRemovalElbowPlot", {"color": "black"}),
    ("doubletRemovalElbowPlot", "doubletRemoval", {"color": "black"}),
    ("doubletRemoval", "addTPsMerge", {"color": "black"}),
    ("addTPsMerge", "SCTransformNormalization", {"color": "black"}),
    ("SCTransformNormalization", "IntegrationDimReduction", {"color": "black"}),
    ("IntegrationDimReduction", "RunUMAP", {"color": "black"}),
    ("RunUMAP", "testDiffClusterResolutions", {"color": "black"}),
    ("RunUMAP", "useChosenClusterResolutions", {"color": "black"}),
    ("useChosenClusterResolutions", "markerDiscovery", {"color": "black"}),
    ("markerDiscovery", "cellCounting", {"color": "black"}),
    ("markerDiscovery", "createShinyApp", {"color": "black"}),
    #nm
    ("metaData", "mt_p1", {"color": "blue"}),
    ("mt_p1", "mt_p2", {"color": "blue"}),
    ("mt_p2", "doubletRemovalElbowPlot", {"color": "blue"}),
    ("doubletRemovalElbowPlot", "doubletRemoval", {"color": "blue"}),
    ("doubletRemoval", "addTPsMerge", {"color": "blue"}),
    ("addTPsMerge", "SCTransformNormalization", {"color": "blue"}),
    ("SCTransformNormalization", "IntegrationDimReduction", {"color": "blue"}),
    ("IntegrationDimReduction", "RunUMAP", {"color": "blue"}),
    ("RunUMAP", "testDiffClusterResolutions", {"color": "blue"}),
    ("RunUMAP", "useChosenClusterResolutions", {"color": "blue"}),
    ("useChosenClusterResolutions", "markerDiscovery", {"color": "blue"}),
    ("markerDiscovery", "cellCounting", {"color": "blue"}),
    ("markerDiscovery", "DGE", {"color": "blue"}),
    ("markerDiscovery", "createShinyApp", {"color": "blue"}),
    #ms
    ("metaData", "mt_p1", {"color": "red"}),
    ("metaData", "multimodalAnalysis", {"color": "red"}),
    ("mt_p1", "mt_p2", {"color": "red"}),
    ("mt_p2", "doubletRemovalElbowPlot", {"color": "red"}),
    ("doubletRemovalElbowPlot", "doubletRemoval", {"color": "red"}),
    ("doubletRemoval", "addTPsMerge", {"color": "red"}),
    ("addTPsMerge", "SCTransformNormalization", {"color": "red"}),
    ("SCTransformNormalization", "IntegrationDimReduction", {"color": "red"}),
    ("IntegrationDimReduction", "RunUMAP", {"color": "red"}),
    ("RunUMAP", "testDiffClusterResolutions", {"color": "red"}),
    ("RunUMAP", "useChosenClusterResolutions", {"color": "red"}),
    ("useChosenClusterResolutions", "multimodalAnalysis", {"color": "red"}),
    ("multimodalAnalysis", "markerDiscovery", {"color": "red"}),
    ("markerDiscovery", "cellCounting", {"color": "red"}),
    ("markerDiscovery", "createShinyApp", {"color": "red"}),
    #mm
    ("metaData", "mt_p1", {"color": "orange"}),
    ("metaData", "multimodalAnalysis", {"color": "orange"}),
    ("mt_p1", "mt_p2", {"color": "orange"}),
    ("mt_p2", "doubletRemovalElbowPlot", {"color": "orange"}),
    ("doubletRemovalElbowPlot", "doubletRemoval", {"color": "orange"}),
    ("doubletRemoval", "addTPsMerge", {"color": "orange"}),
    ("addTPsMerge", "SCTransformNormalization", {"color": "orange"}),
    ("SCTransformNormalization", "IntegrationDimReduction", {"color": "orange"}),
    ("IntegrationDimReduction", "RunUMAP", {"color": "orange"}),
    ("RunUMAP", "testDiffClusterResolutions", {"color": "orange"}),
    ("RunUMAP", "useChosenClusterResolutions", {"color": "orange"}),
    ("useChosenClusterResolutions", "multimodalAnalysis", {"color": "orange"}),
    ("multimodalAnalysis", "markerDiscovery", {"color": "orange"}),
    ("markerDiscovery", "cellCounting", {"color": "orange"}),
    ("markerDiscovery", "DGE", {"color": "orange"}),
    ("markerDiscovery", "createShinyApp", {"color": "orange"}),
    #HTO
    ("demultiplexing", "mt_p1", {"color": "green"}),
    ("mt_p1", "mt_p2", {"color": "green"}),
    ("mt_p2", "doubletRemovalElbowPlot", {"color": "green"}),
    ("doubletRemovalElbowPlot", "addTPsMerge", {"color": "green"}),
    ("addTPsMerge", "SCTransformNormalization", {"color": "green"}),
    ("SCTransformNormalization", "IntegrationDimReduction", {"color": "green"}),
    ("IntegrationDimReduction", "RunUMAP", {"color": "green"}),
    ("RunUMAP", "testDiffClusterResolutions", {"color": "green"}),
    ("RunUMAP", "useChosenClusterResolutions", {"color": "green"}),
    ("useChosenClusterResolutions", "markerDiscovery", {"color": "green"}),
    ("markerDiscovery", "cellCounting", {"color": "green"}),
    ("markerDiscovery", "createShinyApp", {"color": "green"})
])

#pos = nx.get_node_attributes(pipelineRun, "pos")
#color = nx.get_edge_attributes(pipelineRun, "color").values()
#size = list(nx.get_node_attributes(pipelineRun, "size").values())
#print(size)

#nx.draw(pipelineRun, pos, node_size=700, node_shape="h", node_color="w", edge_color=color, with_labels=True, bbox=dict(facecolor="yellow", edgecolor="none", boxstyle="round"))

A = nx.nx_agraph.to_agraph(pipelineRun)
A.layout("dot")
A.draw("test.png")

plt.show()