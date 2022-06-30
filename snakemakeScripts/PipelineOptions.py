#Here you can change certain parameters without fighting through the Snakefile

execTimes = 1

#conda environment
isHHU_HPC = True
def setCondaEnv(onHHU_HPC): #can't access config otherwise here
    global isHHU_HPC 
    isHHU_HPC = onHHU_HPC
if isHHU_HPC:
    countEnv = "envs/HHU_HPC/counting-spec.yml"
    devtoolsEnv = "envs/HHU_HPC/devtools-spec.yml"
    doubletEnv = "envs/HHU_HPC/doublet-spec.yml"
    usualEnv = "envs/HHU_HPC/env-spec.yml"
    markerEnv = "envs/HHU_HPC/marker-spec.yml"
    shinyEnv = "envs/HHU_HPC/shiny-spec.yml"
else:
    #you can insert your customized environments here
    countEnv = "envs/netAccess/counting.yml"
    devtoolsEnv = "envs/netAccess/devtools.yml"
    doubletEnv = "envs/netAccess/doublet.yml"
    usualEnv = "envs/netAccess/env.yml"
    markerEnv = "envs/netAccess/marker.yml"
    shinyEnv = "envs/netAccess/shiny.yml"



#add additional Time and RAM if approximation was not enough
#Time in seconds
addTime = {
    "meta": 0,          #metaData
    "demux": 0,         #demultiplexing
    "mt1": 0,           #mt_p1
    "mt2": 0,           #mt_p2
    "drElbowPlot": 0,   #doubletRemovalElbowPlot
    "doubletRem": 0,    #doubletRemoval
    "addTPs": 0,        #addTPsMerge
    "SCT": 0,           #SCTransformNormalization
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
    "drElbowPlot": 0,   #doubletRemovalElbowPlot
    "doubletRem": 0,    #doubletRemoval
    "addTPs": 0,        #addTPsMerge
    "SCT": 0,           #SCTransformNormalization
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
#Number of cores for rules
#Rules marked with X, don't work with future, more than one cores are useless
numCores = {
    "meta": 1,          # X, metaData
    "demux": 1,         # X, demultiplexing
    "mt1": 1,           # X, mt_p1
    "mt2": 1,           # X, mt_p2
    "drElbowPlot": 1,   # X, doubletRemovalElbowPlot
    "doubletRem": 1,    # X, doubletRemoval, NOT RECOMMENDED, there are issues causing 1 core not to deliver results
    "addTPs": 1,        # X, addTPsMerge
    "SCT": 1,           # X, SCTransformNormalization
    "IntegrDimRed": 1,  #    IntegrationDimReduction
    "UMAP": 1,          # X, RunUMAP
    "testClustRes": 1,  #    testDiffClusterResolutions
    "useClusterRes": 1, # X, useChosenClusterResolutions
    "multimodal": 1,    # X, multimodalAnalysis
    "markerDisc": 1,    #    markerDiscovery
    "cellCount": 1,     # X, cellCounting
    "DGE": 1,           #    DGE
    "shinyApp": 1,      # X, createShinyApp
    "multimodalPlot": 1 # X, multimodalFeaturePlotting
}
