#configfile: "configfiles/config284_215.yaml"

import os
import shutil
import hashlib
import cluster.approximateResources as res
import scripts.SnakefileHelperFuncs as hf


#global variables
doubletEnv = "envs/doublet-spec.yml"# "envs/doublet.yml"
shinyEnv = "envs/shiny-spec.yml" #"envs/shiny.yml"
countEnv = "envs/counting-spec.yml" #"envs/counting.yml"
if(not config["HHU_HPC"]):
  doubletEnv = "envs/doublet.yml"
  shinyEnv = "envs/shiny.yml"
  countEnv = "envs/counting.yml"
projectDirectoryPath = config["projectDirectoryPath"]
if projectDirectoryPath[-1] != "/":
    projectDirectoryPath += "/"
hf.createDirectoriesIfNotExists(projectDirectoryPath)
sampleInputs = hf.transform_sampleInputs(config["sampleInputs"])
num_cells = config["numberOfCells"]
sampleType = "ns"
if config["multiSampled"] and config["multimodal"]:
  sampleType = "mm"
elif config["multiSampled"] and not config["multimodal"]:
  sampleType = "nm"
elif not config["multiSampled"] and config["multimodal"]:
  sampleType = "ms"
additionalTime = 0
assayName = "ADT"
if(config["HTO"]):
  assayName = "HTO"
  sampleType = "HTO"
testClustersName = "plots/" + config["projectName"] + ".res_"


def get_inputs(wildcards):
  inputList = []
  if(hf.checkMinimalInputs(config)):
    cutoff, doublets, conditions = True, True, True
    for i in config["sampleInputs"]:
      if not i["mtCutoff"]:
        cutoff = False
      if not i["dbElbowPlot"] or not i["expectedPercentDoublets"]:
        doublets = False
      if not i["otherMetaData"]:
        conditions = False
    if os.path.isdir(hf.findHash(doubletEnv) + "/lib/R/library/DoubletFinder") and os.path.isdir(hf.findHash(shinyEnv) + "/lib/R/library/ShinyCell") and os.path.isdir(hf.findHash(countEnv) + "/lib/R/library/seurathelpeR"):
      outputStart = projectDirectoryPath + "outputs/" + config["projectName"]
      if config["HTO"]:
        inputList.append(outputStart + ".demux.rds")
      else:
        inputList.append(outputStart + ".meta.rds")
      inputList.append(outputStart + ".mt_p1.rds")
      if cutoff:
        inputList.append(outputStart + ".mt_p2.rds")
        inputList.append(outputStart + ".SCTranDB.rds")
        if doublets:
          if not config["HTO"]:
            inputList.append(outputStart + ".doubR.rds")
          if conditions and config["otherMetaName"]:
            inputList.append(outputStart + ".preprocessedO.rds")
            inputList.append(outputStart + ".normalized.rds")
            inputList.append(outputStart + ".IntDimRed.rds")
            if config["integrationPCs"]:
              inputList.append(outputStart + ".umapped.rds")
              if config["findNeighborsPCs"] and config["choosableResolutions"]:
                if config["chosenResolution"]:
                  inputList.append(outputStart + ".clustered.rds")
                  inputList.append(outputStart + ".markerDisc.rds")
                  if config["multiSampled"]:
                    inputList.append(projectDirectoryPath + "csv/finishedDGE.txt")
                  if config["multimodal"] or config["HTO"]:
                    inputList.append(projectDirectoryPath + "plots/finishedFeatures.txt")
                  inputList.append(projectDirectoryPath + "shinyApp/" + "server.R")
                  inputList.append(projectDirectoryPath + "shinyApp/" + "ui.R")
                  inputList.append(projectDirectoryPath + "shinyApp/howToRunShinyAppOnYourOwnPC.txt")
                  if config["countIdents"]:
                    inputList.append(projectDirectoryPath + "plots/" + config["projectName"] + "." + config["otherMetaName"] + ".barplot.pdf")
                  else:
                    inputList.append(projectDirectoryPath + "workDirectory/countIdentsMissing.txt")
                else:
                  inputList + hf.createMultiSampleInput(projectDirectoryPath, testClustersName, config["choosableResolutions"], ".clusteredDimPlot.pdf")
                  inputList.append(projectDirectoryPath + "workDirectory/chosenResolutionMissing.txt")
              else:
                inputList.append(projectDirectoryPath + "workDirectory/resolutionsNeighPCsMissing.txt")
            else:
              inputList.append(projectDirectoryPath + "workDirectory/intPCsMissing.txt")
          else:
            inputList.append(projectDirectoryPath + "workDirectory/conditionInfosMissing.txt")
        else:
          inputList.append(projectDirectoryPath + "workDirectory/doubletInfosMissing.txt")
      else:
        inputList.append(projectDirectoryPath + "workDirectory/mtCutoffMissing.txt")
    else:
      inputList.append("workDirectory/createDoubletEnv.txt")
      inputList.append("workDirectory/createShinyEnv.txt")
      inputList.append("workDirectory/createCountEnv.txt")
      inputList.append(hf.findHash(doubletEnv) + "/lib/R/library/DoubletFinder")
      inputList.append(hf.findHash(shinyEnv) + "/lib/R/library/ShinyCell")
      inputList.append(hf.findHash(countEnv) + "/lib/R/library/seurathelpeR")
  else:
    print("""Please fill in the minimal amount of inputs needed in the config.yaml:
    'workDirectory', 'rawData', 'projectName', 'projectDirectoryPath', 'multiSampled', 'numberOfCells', 'maxRAM', 'HHU_HPC', 'name' for all samples in 'sampleInputs'
    If you don't know how to fill them in, please check out the 'config_example.yaml'.""")
  return inputList


rule all:
  input:
    get_inputs

rule createDoubletEnv:
  input:
    doubletEnv
  conda:
    doubletEnv
  params:
    mem="2GB",
    time="0:01:59",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="createDoubletEnv"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="createDoubletEnv")
  output:
    temp("workDirectory/createDoubletEnv.txt")
  shell:
    "touch {output}"

rule createCountEnv:
  input:
    countEnv
  conda:
    countEnv
  params:
    mem="2GB",
    time="0:01:59",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="createCountEnv"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="createCountEnv")
  output:
    temp("workDirectory/createCountEnv.txt")
  shell:
    "touch {output}"

rule createShinyEnv:
  input:
    shinyEnv
  conda:
    shinyEnv
  params:
    mem="2GB",
    time="0:01:59",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="createShinyEnv"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="createShinyEnv")
  output:
    temp("workDirectory/createShinyEnv.txt")
  shell:
    "touch {output}"

rule installMissingPackages:
  input:
    expand("workDirectory/createShinyEnv.txt", projectDirPath=projectDirectoryPath),
    expand("workDirectory/createDoubletEnv.txt", projectDirPath=projectDirectoryPath),
    expand("workDirectory/createCountEnv.txt", projectDirPath=projectDirectoryPath)
  params:
    absPath = os.path.realpath("."),
    doublet = hf.findHash(doubletEnv),
    shiny = hf.findHash(shinyEnv),
    countEnv = hf.findHash(countEnv),
    hhuHPC = config["HHU_HPC"],
    mem="2GB",
    time="0:30:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="installMissingPackages"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="installMissingPackages")
  conda:
    "envs/devtools.yml"
  output:
    directory(hf.findHash(doubletEnv) + "/lib/R/library/DoubletFinder"),
    directory(hf.findHash(shinyEnv) + "/lib/R/library/ShinyCell"),
    directory(hf.findHash(countEnv) + "/lib/R/library/seurathelpeR")
  script:
    "scripts/installGithubPacks.R"

rule metaData:
  input:
    expand("{rawData}", rawData=config["rawData"])
  params:
    projectDirPath=projectDirectoryPath,
    names = sampleInputs["name"],
    project = config["projectName"],
    time=res.approxWalltime("meta", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("meta", sampleType, num_cells, 2),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="metaData"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="metaData")
  conda:
    "envs/env.yml"
  output:
    expand(["{projectDirPath}outputs/{project}.meta.rds", "{projectDirPath}outputs/{project}.rawData.rds"], project=config["projectName"], projectDirPath=projectDirectoryPath) if config["multimodal"] else expand("{projectDirPath}outputs/{project}.meta.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
    #hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleInputs["name"], ".meta.rds")
  script:
    "scripts/MetaData.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/1_meta.txt", 3)

rule demultiplexing:
  input:
    expand("{rawData}", rawData=config["rawData"])
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    assayname = assayName,
    names = sampleInputs["name"],
    time=res.approxWalltime("demultiplex", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("demultiplex", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="demultiplexing"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="demultiplexing")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}outputs/{project}.demux.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
    #hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleInputs["name"], ".demux.rds")
  script:
    "scripts/demultiplexing.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/97_demux_umap.txt", 3)

rule mt_p1:
  input:
    expand("{projectDirPath}outputs/{project}.meta.rds", project=config["projectName"], projectDirPath=projectDirectoryPath) if not config["HTO"] else expand("{projectDirPath}outputs/{project}.demux.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    names = sampleInputs["name"],
    pattern = config["pattern"],
    time=res.approxWalltime("mt1", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("mt1", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="mt_p1"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="mt_p1")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}outputs/{project}.mt_p1.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/mt_p1.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/2_mt1.txt", 3)

rule mt_p2:
  input: 
    expand("{projectDirPath}outputs/{project}.mt_p1.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    samples = sampleInputs["mtCutoff"],
    names = sampleInputs["name"],
    time=res.approxWalltime("mt2", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("mt2", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="mt_p2"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="mt_p2")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}outputs/{project}.mt_p2.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/mt_p2.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/3_mt2.txt", 3)

rule doubletRemovalElbowPlot:
  input:
    expand("{projectDirPath}outputs/{project}.mt_p2.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    names = sampleInputs["name"],
    time=res.approxWalltime("dbElb", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("dbElb", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="doubletRemovalElbowPlot"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="doubletRemovalElbowPlot")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}outputs/{project}.SCTranDB.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/DBElbowPlotter.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/4_dbElbow.txt", 3)

rule doubletRemoval:
  input:
    expand("{projectDirPath}outputs/{project}.SCTranDB.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    names = sampleInputs["name"],
    roundingValues = sampleInputs["expectedPercentDoublets"],
    dim_PCs = sampleInputs["dbElbowPlot"],
    time=res.approxWalltime("dbRem", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("dbRem", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="doubletRemoval"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="doubletRemoval")
  conda:
    doubletEnv
  output:
    expand("{projectDirPath}outputs/{project}.doubR.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/DoubletRemoval.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/5_dbRem.txt", 3)

rule addTPsMerge:
  input:
    expand("{projectDirPath}outputs/{project}.SCTranDB.rds", project=config["projectName"], projectDirPath=projectDirectoryPath) if config["HTO"] else expand("{projectDirPath}outputs/{project}.doubR.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
    #hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleInputs["name"], ".doubR.rds")
  params:
    projectDirPath = projectDirectoryPath,
    meta = sampleInputs["otherMetaData"],
    metaName = config["otherMetaName"],
    time=res.approxWalltime("merge", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("merge", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="addTPsMerge"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="addTPsMerge")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}outputs/{project}.preprocessedO.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
    #hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleInputs["name"], ".preprocessed.rds")
  script:
    "scripts/addMetaAndMerge.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/6_merge.txt", 3)
  
rule SCTransformNormalization:
  input:
    expand("{projectDirPath}outputs/{project}.preprocessedO.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    time=res.approxWalltime("sct", sampleType, num_cells),
    mem=res.approxRAM("sct", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="SCTransformNormalization"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="SCTransformNormalization")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}outputs/{project}.normalized.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/SCTraNormalisation.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/7_sctran.txt", 3)

rule IntegrationDimReduction:
  input:
    expand("{projectDirPath}outputs/{project}.normalized.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
    #hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleInputs["name"], ".normalized.rds")
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    maxRAM = config["maxRAM"],
    time=res.approxWalltime("integration", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("integration", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="IntegrationDimReduction"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="IntegrationDimReduction")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}outputs/{project}.IntDimRed.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/IntegrationDimReduction.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/8_integration.txt", 3)

rule RunUMAP:
  input:
    expand("{projectDirPath}outputs/{project}.IntDimRed.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    integrationPCs = config["integrationPCs"],
    time=res.approxWalltime("UMAP", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("UMAP", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="RunUMAP"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="RunUMAP")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}outputs/{project}.umapped.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/RunUMAP.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/9_UMAP.txt", 3)

rule testDiffClusterResolutions:
  input:
    expand("{projectDirPath}outputs/{project}.umapped.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    findNeighborsPCs = config["findNeighborsPCs"],
    resolutions = config["choosableResolutions"],
    time=res.approxWalltime("tCluster", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("tCluster", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="testDiffClusterResolutions"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="testDiffClusterResolutions")
  conda:
    "envs/env.yml"
  output:
    #expand("{projectDirPath}plots/finished.txt", project=config["projectName"], projectDirPath=projectDirectoryPath)
    hf.createMultiSampleInput(projectDirectoryPath, testClustersName, config["choosableResolutions"], ".clusteredDimPlot.pdf")
  script:
    "scripts/testDiffClusterRes.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/9_3_testcluster.txt", 3)

rule useChosenClusterResolutions:
  input:
    expand("{projectDirPath}outputs/{project}.umapped.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    findNeighborsPCs = config["findNeighborsPCs"],
    resolution = config["chosenResolution"],
    condition = config["otherMetaName"],
    time=res.approxWalltime("cCluster", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("cCluster", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="useChosenClusterResolutions"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="useChosenClusterResolutions")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}outputs/{project}.clustered.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/useChosenClusterRes.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/90_chosenCluster.txt", 3)

rule multimodalAnalysis:
  input:
    expand(["{projectDirPath}outputs/{project}.clustered.rds", "{projectDirPath}outputs/{project}.rawData.rds"], project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    assayname = assayName,
    time=res.approxWalltime("mmodalAssay", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("mmodalAssay", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="multimodalAnalysis"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="multimodalAnalysis")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}outputs/{project}.multimodal.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/addAssayData.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/81_addAssay.txt", 3)

rule markerDiscovery:
  input:
    expand("{projectDirPath}outputs/{project}.multimodal.rds", project=config["projectName"], projectDirPath=projectDirectoryPath) if config["multimodal"] and not config["HTO"] else expand("{projectDirPath}outputs/{project}.clustered.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
    #expand(["projects/{projectDirName}outputs/{project}.clustered.rds", "projects/{projectDirName}outputs/{projects}.rawData.rds"], project=config["projectName"], projectDirName=config["projectDirName"]) if config["mutlimodal"] else "projects/{projectDirName}outputs/{project}.clustered.rds"
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    time=res.approxWalltime("marker", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("marker", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="markerDiscovery"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="markerDiscovery")
  conda:
    "envs/marker.yml"
  output:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/markerDiscovery.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/91_marker.txt", 3)

rule cellCounting:
  input:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    countIdents = config["otherMetaName"],
    resolution = config["chosenResolution"],
    multiSampled = config["multiSampled"],
    project = config["projectName"],
    time=res.approxWalltime("count", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("count", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="cellCounting"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="cellCounting")
  conda:
    countEnv
  output:
    expand("{projectDirPath}plots/{project}.{condition}.barplot.pdf", project=config["projectName"], projectDirPath=projectDirectoryPath, condition=config["otherMetaName"])
  script:
    "scripts/cellCounting.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/92_ncCount.txt", 3)

rule DGE:
  input:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    metaCondition = config["otherMetaName"],
    metaIdent = sampleInputs["otherMetaData"],
    resolution = config["chosenResolution"],
    time=res.approxWalltime("DGE", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("DGE", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="DGE"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="DGE")
  conda:
    "envs/marker.yml"
  output:
    expand("{projectDirPath}csv/finishedDGE.txt", projectDirPath=projectDirectoryPath)
  script:
    "scripts/DGE.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/94_dge.txt", 3)

rule createShinyApp:
  input:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    multiSampled = config["multiSampled"],
    time=res.approxWalltime("shiny", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("shiny", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="createShinyApp"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="createShinyApp")
  conda:
    shinyEnv
  output:
    expand("{projectDirPath}shinyApp/server.R", project=config["projectName"], projectDirPath=projectDirectoryPath),
    expand("{projectDirPath}shinyApp/ui.R", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/createShinyApp.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/95_shiny.txt", 3)

rule multimodalFeaturePlotting:
  input:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    assayname = assayName,
    time=res.approxWalltime("mmodalplot", sampleType, num_cells, additionalTime),
    mem=res.approxRAM("mmodalplot", sampleType, num_cells),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="multimodalFeaturePlotting"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="multimodalFeaturePlotting")
  conda:
    "envs/env.yml"
  output:
    expand("{projectDirPath}plots/finishedFeatures.txt", project=config["projectName"], projectDirPath=projectDirectoryPath)
  script:
    "scripts/multimodalFeaturePlotting.R"
#  benchmark:
#    repeat("benchmarks/benchmark284_215_multimod/96_multimodalplot.txt", 3)

rule copyShinyAppInstructions:
  params:
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="copyInstructions"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="copyInstructions")
  output:
    expand("{projectDirPath}shinyApp/howToRunShinyAppOnYourOwnPC.txt", projectDirPath=projectDirectoryPath)
  run:
    shutil.copyfile("workDirectory/howToRunShinyAppOnYourOwnPC.txt", projectDirectoryPath + "shinyApp/howToRunShinyAppOnYourOwnPC.txt")

####################################  missing paramter rules  ###############################################


rule missingMTCutoff:
  input:
    expand("{projectDirPath}outputs/{project}.mt_p1.rds", project=config["projectName"], projectDirPath=projectDirectoryPath),
    "errorMessage/mtCutoffMissing.txt"
  params:
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingMTCutoff"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingMTCutoff")
  output:
    temp(expand("{projectDirPath}workDirectory/mtCutoffMissing.txt", projectDirPath=projectDirectoryPath))
  run:
    hf.raiseInputMissingException("errorMessage/mtCutoffMissing.txt")
    shell("touch {output}")

rule missingDoubletInfos:
  input:
    expand("{projectDirPath}outputs/{project}.SCTranDB.rds", project=config["projectName"], projectDirPath=projectDirectoryPath),
    "errorMessage/doubletInfosMissing.txt"
  params:
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingDoubletInfos"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingDoubletInfos")
  output:
    temp(expand("{projectDirPath}workDirectory/doubletInfosMissing.txt", projectDirPath=projectDirectoryPath))
  run:
    hf.raiseInputMissingException("errorMessage/doubletInfosMissing.txt")
    shell("touch {output}")

rule missingConditions:
  input:
    expand("{projectDirPath}outputs/{project}.doubR.rds", project=config["projectName"], projectDirPath=projectDirectoryPath),
    "errorMessage/conditionInfosMissing.txt"
  params:
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingConditions"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingConditions")
  output:
    temp(expand("{projectDirPath}workDirectory/conditionInfosMissing.txt", projectDirPath=projectDirectoryPath))
  run:
    hf.raiseInputMissingException("errorMessage/conditionInfosMissing.txt")
    shell("touch {output}")

rule missingIntPCs:
  input:
    expand("{projectDirPath}outputs/{project}.IntDimRed.rds", project=config["projectName"], projectDirPath=projectDirectoryPath),
    "errorMessage/intPCsMissing.txt"
  params:
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingIntPCs"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingIntPCs")
  output:
    temp(expand("{projectDirPath}workDirectory/intPCsMissing.txt", projectDirPath=projectDirectoryPath))
  run:
    hf.raiseInputMissingException("errorMessage/intPCsMissing.txt")
    shell("touch {output}")

rule missingResolutionsNeighPC:
  input:
    expand("{projectDirPath}outputs/{project}.umapped.rds", project=config["projectName"], projectDirPath=projectDirectoryPath),
    "errorMessage/resolutionsNeighPCsMissing.txt"
  params:
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingResolutionsNeighPC"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingResolutionsNeighPC")
  output:
    temp(expand("{projectDirPath}workDirectory/resolutionsNeighPCsMissing.txt", projectDirPath=projectDirectoryPath))
  run:
    hf.raiseInputMissingException("errorMessage/resolutionsNeighPCsMissing.txt")
    shell("touch {output}")

rule missingChosenResolution:
  input:
    hf.createMultiSampleInput(projectDirectoryPath, testClustersName, config["choosableResolutions"], ".clusteredDimPlot.pdf"),
    "errorMessage/chosenResolutionMissing.txt"
  params:
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingChosenResolution"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingChosenResolution")
  output:
    temp(expand("{projectDirPath}workDirectory/chosenResolutionMissing.txt", projectDirPath=projectDirectoryPath))
  run:
    hf.raiseInputMissingException("errorMessage/chosenResolutionMissing.txt")
    shell("touch {output}")

rule missingCountIdents:
  input:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath),
    "errorMessage/countIdentsMissing.txt"
  params:
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingCountIdents"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingCountIdents")
  output:
    temp(expand("{projectDirPath}workDirectory/countIdentsMissing.txt", projectDirPath=projectDirectoryPath))
  run:
    hf.raiseInputMissingException("errorMessage/countIdentsMissing.txt")
    shell("touch {output}")
