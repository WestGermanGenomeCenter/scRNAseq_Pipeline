#Here you can change certain parameters without fighting through the Snakefile

#conda environment
countEnv = "envs/counting-spec.yml"     # "envs/counting"
devtoolsEnv = "envs/devtools-spec.yml"  # "envs/devtools"
doubletEnv = "envs/doublet-spec.yml"    # "envs/doublet"
usualEnv = "envs/env-spec.yml"          # "envs/env"
markerEnv = "envs/marker-spec.yml"      # "envs/marker"
shinyEnv = "envs/shiny-spec.yml"        # "envs/shiny"

#add additional Time and RAM if approximation was not enough
#Time in seconds
addTime = {
    "meta": 0,          #metaData
    "demux": 0,         #demultiplexing
    "mt1": 0,           #mt_p1
    "mt2": 0,           #mt_p2
    "drElbowPlot": 0,   #doubletRemovalElbowPlot
    "doubletRem": 0,    #doubletRemoval
    "addTPs": 600,        #addTPsMerge
    "SCT": 300,           #SCTransformNormalization
    "IntegrDimRed": 0,  #IntegrationDimReduction
    "UMAP": 0,          #RunUMAP
    "testClustRes": 0,  #testDiffClusterResolutions
    "useClusterRes": 0, #useChosenClusterResolutions
    "multimodal": 0,    #multimodalAnalysis
    "markerDisc": 0,    #markerDiscovery
    "cellCount": 0,     #cellCounting
    "DGE": 0,           #DGE
    "shinyApp": 0,      #createShinyApp
    "multimodalPlot": 0 #multimodalFeaturePlotting
}
#RAM in GB
addRAM = {
    "meta": 0,          #metaData
    "demux": 0,         #demultiplexing
    "mt1": 0,           #mt_p1
    "mt2": 0,           #mt_p2
    "drElbowPlot": 20,   #doubletRemovalElbowPlot
    "doubletRem": 0,    #doubletRemoval
    "addTPs": 10,        #addTPsMerge
    "SCT": 10,           #SCTransformNormalization
    "IntegrDimRed": 0,  #IntegrationDimReduction
    "UMAP": 0,          #RunUMAP
    "testClustRes": 0,  #testDiffClusterResolutions
    "useClusterRes": 0, #useChosenClusterResolutions
    "multimodal": 0,    #multimodalAnalysis
    "markerDisc": 0,    #markerDiscovery
    "cellCount": 0,     #cellCounting
    "DGE": 0,           #DGE
    "shinyApp": 0,      #createShinyApp
    "multimodalPlot": 0 #multimodalFeaturePlotting
}