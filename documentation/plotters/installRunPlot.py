import networkx as nx

installRun = nx.MultiDiGraph()
installRun.add_nodes_from([
    ("DoubletFinder", {"fontname": "calibri", "shape": "folder", "style": "filled"}),
    ("ShinyCell", {"fontname": "calibri", "shape": "folder", "style": "filled"}),
    ("seurathelpeR", {"fontname": "calibri", "shape": "folder", "style": "filled"}),
    ("GenomeInfoDb", {"fontname": "calibri", "shape": "folder", "style": "dashed"}),
    ("GenomeInfoDbData", {"fontname": "calibri", "shape": "folder", "style": "dashed"}),
    ("GenomicRanges", {"fontname": "calibri", "shape": "folder", "style": "dashed"}),
    ("SummarizedExperiment", {"fontname": "calibri", "shape": "folder", "style": "dashed"}),
    ("SingleCellExperiment", {"fontname": "calibri", "shape": "folder", "style": "dashed"}),
    ("createDoubletEnv", {"fontname": "calibri"}),
    ("createShinyEnv", {"fontname": "calibri"}),
    ("createCountEnv", {"fontname": "calibri"}),
    ("installMissingPackages", {"fontname": "calibri"}),
])
installRun.add_edges_from([
    ("DoubletFinder", "createDoubletEnv", {"color": "black"}),
    ("ShinyCell", "createShinyEnv", {"color": "black"}),
    ("GenomeInfoDb", "createCountEnv", {"color": "black"}),
    ("GenomeInfoDbData", "createCountEnv", {"color": "black"}),
    ("GenomicRanges", "createCountEnv", {"color": "black"}),
    ("SummarizedExperiment", "createCountEnv", {"color": "black"}),
    ("SingleCellExperiment", "createCountEnv", {"color": "black"}),
    ("seurathelpeR", "createCountEnv", {"color": "black"}),
    ("createDoubletEnv", "installMissingPackages", {"color": "black"}),
    ("createShinyEnv", "installMissingPackages", {"color": "black"}),
    ("createCountEnv", "installMissingPackages", {"color": "black"}),
])
installRun.add_nodes_from([
    ("legend:", {"fontname": "calibri", "shape": "none"}),
    ("package\nnot on conda", {"fontname": "calibri", "shape": "folder", "style": "filled"}),
    ("issues with\nHPC's conda", {"fontname": "calibri", "shape": "folder", "style": "dashed"}),
])
installRun.add_edges_from([
    ("legend:", "package\nnot on conda", {"color": "white"}),
    ("legend:", "issues with\nHPC's conda", {"color": "white"}),
])
installRun.graph["graph"] = {"rankdir": "LR",
                             "size": "8.00,10.00"}
