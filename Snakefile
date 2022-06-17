#configfile: "configfiles/config284_215.yaml"

import os
import shutil
import hashlib
import snakemakeScripts.cluster.approximateResources as res
import snakemakeScripts.SnakefileHelperFuncs as hf
import snakemakeScripts.PipelineOptions as pOpt


projectDirectoryPath = config["projectDirectoryPath"]
if projectDirectoryPath[-1] != "/":
    projectDirectoryPath += "/"
hf.createDirectoriesIfNotExists(projectDirectoryPath)
sampleInputs, otherMetaData = hf.transform_sampleInputs(config["sampleInputs"])
sampleNames = list(sampleInputs.keys())
num_cells = config["numberOfCells"]
sampleType = "ns"
if config["multiSampled"] and config["multimodal"]:
  sampleType = "mm"
elif config["multiSampled"] and not config["multimodal"]:
  sampleType = "nm"
elif not config["multiSampled"] and config["multimodal"]:
  sampleType = "ms"
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
    if os.path.isdir(hf.findHash(pOpt.doubletEnv) + "/lib/R/library/DoubletFinder") and os.path.isdir(hf.findHash(pOpt.shinyEnv) + "/lib/R/library/ShinyCell") and os.path.isdir(hf.findHash(pOpt.countEnv) + "/lib/R/library/seurathelpeR"):
      outputStart = projectDirectoryPath + "outputs/" + config["projectName"]
      if config["HTO"]:
        inputList = hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".demux.rds")
      else:
        if config["multimodal"]:
          inputList = hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".meta.rds", config["projectName"])
        else:
          inputList = hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".meta.rds")
      inputList = inputList + hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".mt_p1.rds")
      if cutoff:
        inputList = inputList + hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".mt_p2.rds")
        inputList = inputList + hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".SCTranDB.rds")
        if doublets:
          if not config["HTO"]:
            inputList = inputList + hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".doubR.rds")
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
                    for i in hf.createCombinations(config["otherMetaName"], combiOnly=False):
                      inputList.append(projectDirectoryPath + "csv/" + i + "/" + i + ".finishedDGE.txt")
                  if config["multimodal"] or config["HTO"]:
                    inputList.append(projectDirectoryPath + "plots/finishedFeatures.txt")
                  inputList.append(projectDirectoryPath + "shinyApp/" + "server.R")
                  inputList.append(projectDirectoryPath + "shinyApp/" + "ui.R")
                  inputList.append(projectDirectoryPath + "shinyApp/howToRunShinyAppOnYourOwnPC.txt")
                  inputList += hf.createMultiMetaCountInput(projectDirectoryPath, "plots/" + config["projectName"] + ".", config["otherMetaName"], ".barplot.pdf")
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
      inputList.append(hf.findHash(pOpt.doubletEnv) + "/lib/R/library/DoubletFinder")
      inputList.append(hf.findHash(pOpt.shinyEnv) + "/lib/R/library/ShinyCell")
      inputList.append(hf.findHash(pOpt.countEnv) + "/lib/R/library/seurathelpeR")
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
    pOpt.doubletEnv
  conda:
    pOpt.doubletEnv
  params:
    cores=1,
    mem="2GB",
    time="0:01:59",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="createDoubletEnv"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="createDoubletEnv")
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="createDoubletEnv")
  output:
    temp("workDirectory/createDoubletEnv.txt")
  shell:
    "touch {output}"


rule createCountEnv:
  input:
    pOpt.countEnv
  conda:
    pOpt.countEnv
  params:
    cores=1,
    mem="2GB",
    time="0:01:59",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="createCountEnv"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="createCountEnv")
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="createCountEnv")
  output:
    temp("workDirectory/createCountEnv.txt")
  shell:
    "touch {output}"


rule createShinyEnv:
  input:
    pOpt.shinyEnv
  conda:
    pOpt.shinyEnv
  params:
    cores=1,
    mem="2GB",
    time="0:01:59",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="createShinyEnv"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="createShinyEnv")
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="createShinyEnv")
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
    doublet = hf.findHash(pOpt.doubletEnv),
    shiny = hf.findHash(pOpt.shinyEnv),
    count = hf.findHash(pOpt.countEnv),
    hhuHPC = config["HHU_HPC"],
    cores=1,
    mem="2GB",
    time="0:30:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="installMissingPackages"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="installMissingPackages")
  conda:
    pOpt.devtoolsEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="installMissingPackages")
  output:
    directory(hf.findHash(pOpt.doubletEnv) + "/lib/R/library/DoubletFinder"),
    directory(hf.findHash(pOpt.shinyEnv) + "/lib/R/library/ShinyCell"),
    directory(hf.findHash(pOpt.countEnv) + "/lib/R/library/seurathelpeR")
  script:
    "scripts/installGithubPacks.R"


rule metaData:
  input:
    expand("{rawData}", rawData=config["rawData"])
  params:
    projectDirPath=projectDirectoryPath,
    names = sampleNames,
    project = config["projectName"],
    cores=pOpt.numCores["meta"],
    time=res.approxWalltime("meta", sampleType, num_cells, additionalTime=pOpt.addTime["meta"]),
    mem=res.approxRAM("meta", sampleType, num_cells, additionalRAM=pOpt.addRAM["meta"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="metaData"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="metaData")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="metaData")
  output:
    hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".meta.rds", project=config["projectName"]) if config["multimodal"] else hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".meta.rds")
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/00_meta.txt", 3)
  script:
    "scripts/MetaData.R"


rule demultiplexing:
  input:
    expand("{rawData}", rawData=config["rawData"])
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    assayname = assayName,
    names = sampleNames,
    cores=pOpt.numCores["demux"],
    time=res.approxWalltime("demultiplex", sampleType, num_cells, additionalTime=pOpt.addTime["demux"]),
    mem=res.approxRAM("demultiplex", sampleType, num_cells, additionalRAM=pOpt.addRAM["demux"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="demultiplexing"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="demultiplexing")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="demultiplexing")
  output:
    hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".demux.rds")
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/01_demux.txt", 3)
  script:
    "scripts/demultiplexing.R"


rule mt_p1:
  input:
    expand("{projectDirPath}outputs/{{names}}.meta.rds", projectDirPath=projectDirectoryPath) if not config["HTO"] else expand("{projectDirPath}outputs/{{names}}.demux.rds", projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    names ="{names}",
    pattern = config["pattern"],
    cores=pOpt.numCores["mt1"],
    time=res.approxWalltime("mt1", sampleType, num_cells, additionalTime=pOpt.addTime["mt1"]),
    mem=res.approxRAM("mt1", sampleType, num_cells, additionalRAM=pOpt.addRAM["mt1"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.{{names}}.errors", projectDirPath=projectDirectoryPath, rule="mt_p1"),
    output=expand("{projectDirPath}clusterLogs/{rule}.{{names}}.output", projectDirPath=projectDirectoryPath, rule="mt_p1")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.{{names}}.log", projectDirPath=projectDirectoryPath, rule="mt_p1")
  output:
    expand("{projectDirPath}outputs/{{names}}.mt_p1.rds", projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/02_mt1.{{names}}.txt", 3)
  script:
    "scripts/mt_p1.R"


rule mt_p2:
  input: 
    expand("{projectDirPath}outputs/{{names}}.mt_p1.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    names ="{names}",
    samples = lambda wcs: sampleInputs[wcs.names]['mtCutoff'],
    cores=pOpt.numCores["mt2"],
    time=res.approxWalltime("mt2", sampleType, num_cells, additionalTime=pOpt.addTime["mt2"]),
    mem=res.approxRAM("mt2", sampleType, num_cells, additionalRAM=pOpt.addRAM["mt2"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.{{names}}.errors", projectDirPath=projectDirectoryPath, rule="mt_p2"),
    output=expand("{projectDirPath}clusterLogs/{rule}.{{names}}.output", projectDirPath=projectDirectoryPath, rule="mt_p2")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.{{names}}.log", projectDirPath=projectDirectoryPath, rule="mt_p2")
  output:
    expand("{projectDirPath}outputs/{{names}}.mt_p2.rds", projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/03_mt2.{{names}}.txt", 3)
  script:
    "scripts/mt_p2.R"


rule doubletRemovalElbowPlot:
  input:
    expand("{projectDirPath}outputs/{{names}}.mt_p2.rds", projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    names = "{names}",
    cores=pOpt.numCores["drElbowPlot"],
    time=res.approxWalltime("dbElb", sampleType, num_cells, additionalTime=pOpt.addTime["drElbowPlot"]),
    mem=res.approxRAM("dbElb", sampleType, num_cells, additionalRAM=pOpt.addRAM["drElbowPlot"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.{{names}}.errors", projectDirPath=projectDirectoryPath, rule="doubletRemovalElbowPlot"),
    output=expand("{projectDirPath}clusterLogs/{rule}.{{names}}.output", projectDirPath=projectDirectoryPath, rule="doubletRemovalElbowPlot")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.{{names}}.log", projectDirPath=projectDirectoryPath, rule="doubletRemovalElbowPlot")
  output:
    expand("{projectDirPath}outputs/{{names}}.SCTranDB.rds", projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/04_dbElb.{{names}}.txt", 3)
  script:
    "scripts/DBElbowPlotter.R"


rule doubletRemoval:
  input:
    expand("{projectDirPath}outputs/{{names}}.SCTranDB.rds", projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    names = "{names}",
    roundingValues = lambda wcs: sampleInputs[wcs.names]['expectedPercentDoublets'],
    dim_PCs = lambda wcs: sampleInputs[wcs.names]['dbElbowPlot'],
    cores=pOpt.numCores["doubletRem"],
    time=res.approxWalltime("dbRem", sampleType, num_cells, additionalTime=pOpt.addTime["doubletRem"]),
    mem=res.approxRAM("dbRem", sampleType, num_cells, additionalRAM=pOpt.addRAM["doubletRem"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.{{names}}.errors", projectDirPath=projectDirectoryPath, rule="doubletRemoval"),
    output=expand("{projectDirPath}clusterLogs/{rule}.{{names}}.output", projectDirPath=projectDirectoryPath, rule="doubletRemoval")
  conda:
    pOpt.doubletEnv
  log:
    expand("{projectDirPath}logs/{rule}.{{names}}.log", projectDirPath=projectDirectoryPath, rule="doubletRemoval")
  output:
    expand("{projectDirPath}outputs/{{names}}.doubR.rds", projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/05_dbRem.{{names}}.txt", 3)
  script:
    "scripts/DoubletRemoval.R"


rule addTPsMerge:
  input:
    hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".SCTranDB.rds") if config["HTO"] else hf.createMultiSampleInput(projectDirectoryPath, "outputs/", sampleNames, ".doubR.rds")
  params:
    projectDirPath = projectDirectoryPath,
    meta = otherMetaData,
    metaName = config["otherMetaName"],
    condCombi = hf.createCombinations(config["otherMetaName"]),
    cores=pOpt.numCores["addTPs"],
    time=res.approxWalltime("merge", sampleType, num_cells, additionalTime=pOpt.addTime["addTPs"]),
    mem=res.approxRAM("merge", sampleType, num_cells, additionalRAM=pOpt.addRAM["addTPs"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="addTPsMerge"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="addTPsMerge")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="addTPsMerge")
  output:
    expand("{projectDirPath}outputs/{project}.preprocessedO.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/06_merge.txt", 3)
  script:
    "scripts/addMetaAndMerge.R"
  

rule SCTransformNormalization:
  input:
    expand("{projectDirPath}outputs/{project}.preprocessedO.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    cores=pOpt.numCores["SCT"],
    time=res.approxWalltime("sct", sampleType, num_cells, additionalTime=pOpt.addTime["SCT"]),
    mem=res.approxRAM("sct", sampleType, num_cells, additionalRAM=pOpt.addRAM["SCT"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="SCTransformNormalization"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="SCTransformNormalization")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="SCTransformNormalization")
  output:
    expand("{projectDirPath}outputs/{project}.normalized.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/07_SCT.txt", 3)
  script:
    "scripts/SCTraNormalisation.R"


rule IntegrationDimReduction:
  input:
    expand("{projectDirPath}outputs/{project}.normalized.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    maxRAM = config["maxRAM"],
    cores = pOpt.numCores["IntegrDimRed"],
    time=res.approxWalltime("integration", sampleType, num_cells, additionalTime=pOpt.addTime["IntegrDimRed"]),
    mem=res.approxRAM("integration", sampleType, num_cells, additionalRAM=pOpt.addRAM["IntegrDimRed"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="IntegrationDimReduction"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="IntegrationDimReduction")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="IntegrationDimReduction")
  output:
    expand("{projectDirPath}outputs/{project}.IntDimRed.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/08_integDimRed.txt", 3)
  script:
    "scripts/IntegrationDimReduction.R"


rule RunUMAP:
  input:
    expand("{projectDirPath}outputs/{project}.IntDimRed.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    integrationPCs = config["integrationPCs"],
    cores = pOpt.numCores["UMAP"],
    time=res.approxWalltime("UMAP", sampleType, num_cells, additionalTime=pOpt.addTime["UMAP"]),
    mem=res.approxRAM("UMAP", sampleType, num_cells, additionalRAM=pOpt.addRAM["UMAP"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="RunUMAP"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="RunUMAP")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="RunUMAP")
  output:
    expand("{projectDirPath}outputs/{project}.umapped.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/09_umap.txt", 3)
  script:
    "scripts/RunUMAP.R"


rule testDiffClusterResolutions:
  input:
    expand("{projectDirPath}outputs/{project}.umapped.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    findNeighborsPCs = config["findNeighborsPCs"],
    resolutions = config["choosableResolutions"],
    cores = pOpt.numCores["testClustRes"],
    isMultiSampled = config["multiSampled"],
    time=res.approxWalltime("tCluster", sampleType, num_cells, additionalTime=pOpt.addTime["testClustRes"]),
    mem=res.approxRAM("tCluster", sampleType, num_cells, additionalRAM=pOpt.addRAM["testClustRes"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="testDiffClusterResolutions"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="testDiffClusterResolutions")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="testDiffClusterResolutions")
  output:
    hf.createMultiSampleInput(projectDirectoryPath, testClustersName, config["choosableResolutions"], ".clusteredDimPlot.pdf")
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/10_testRes.txt", 3)
  script:
    "scripts/testDiffClusterRes.R"


rule useChosenClusterResolutions:
  input:
    expand("{projectDirPath}outputs/{project}.umapped.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    findNeighborsPCs = config["findNeighborsPCs"],
    resolution = config["chosenResolution"],
    condition = config["otherMetaName"],
    cores = pOpt.numCores["useClusterRes"],
    time=res.approxWalltime("cCluster", sampleType, num_cells, additionalTime=pOpt.addTime["useClusterRes"]),
    mem=res.approxRAM("cCluster", sampleType, num_cells, additionalRAM=pOpt.addRAM["useClusterRes"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="useChosenClusterResolutions"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="useChosenClusterResolutions")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="useChosenClusterResolutions")
  output:
    expand("{projectDirPath}outputs/{project}.clustered.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/11_useRes.txt", 3)
  script:
    "scripts/useChosenClusterRes.R"


rule multimodalAnalysis:
  input:
    expand(["{projectDirPath}outputs/{project}.clustered.rds", "{projectDirPath}outputs/{project}.rawData.rds"], project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    assayname = assayName,
    cores = pOpt.numCores["multimodal"],
    time=res.approxWalltime("mmodalAssay", sampleType, num_cells, additionalTime=pOpt.addTime["multimodal"]),
    mem=res.approxRAM("mmodalAssay", sampleType, num_cells, additionalRAM=pOpt.addRAM["multimodal"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="multimodalAnalysis"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="multimodalAnalysis")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="multimodalAnalysis")
  output:
    expand("{projectDirPath}outputs/{project}.multimodal.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/12_addAssay.txt", 3)
  script:
    "scripts/addAssayData.R"


rule markerDiscovery:
  input:
    expand("{projectDirPath}outputs/{project}.multimodal.rds", project=config["projectName"], projectDirPath=projectDirectoryPath) if config["multimodal"] and not config["HTO"] else expand("{projectDirPath}outputs/{project}.clustered.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    cores = pOpt.numCores["markerDisc"],
    time=res.approxWalltime("marker", sampleType, num_cells, additionalTime=pOpt.addTime["markerDisc"]),
    mem=res.approxRAM("marker", sampleType, num_cells, additionalRAM=pOpt.addRAM["markerDisc"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="markerDiscovery"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="markerDiscovery")
  conda:
    pOpt.markerEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="markerDiscovery")
  output:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/13_markerDisc.txt", 3)
  script:
    "scripts/markerDiscovery.R"


rule cellCounting:
  input:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    countCondition = "{condition}",
    resolution = config["chosenResolution"],
    multiSampled = config["multiSampled"],
    project = config["projectName"],
    cores = pOpt.numCores["cellCount"],
    time=res.approxWalltime("count", sampleType, num_cells, additionalTime=pOpt.addTime["cellCount"]),
    mem=res.approxRAM("count", sampleType, num_cells, additionalRAM=pOpt.addRAM["cellCount"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.{{condition}}.errors", projectDirPath=projectDirectoryPath, rule="cellCounting"),
    output=expand("{projectDirPath}clusterLogs/{rule}.{{condition}}.output", projectDirPath=projectDirectoryPath, rule="cellCounting")
  conda:
    pOpt.countEnv
  log:
    expand("{projectDirPath}logs/{rule}.{{condition}}.log", projectDirPath=projectDirectoryPath, rule="cellCounting")
  output:
    expand("{projectDirPath}plots/{project}.{{condition}}.barplot.pdf", project=config["projectName"], projectDirPath=projectDirectoryPath, condition=config["otherMetaName"])
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/14_cCount.{{condition}}.txt", 3)
  script:
    "scripts/cellCounting.R"


rule DGE:
  input:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    metaCondition = "{condition}",
    resolution = config["chosenResolution"],
    cores = pOpt.numCores["DGE"],
    time=res.approxWalltime("DGE", sampleType, num_cells, additionalTime=pOpt.addTime["DGE"]),
    mem=res.approxRAM("DGE", sampleType, num_cells, additionalRAM=pOpt.addRAM["DGE"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.{{condition}}.errors", projectDirPath=projectDirectoryPath, rule="DGE"),
    output=expand("{projectDirPath}clusterLogs/{rule}.{{condition}}.output", projectDirPath=projectDirectoryPath, rule="DGE")
  conda:
    pOpt.markerEnv
  log:
    expand("{projectDirPath}logs/{rule}.{{condition}}.log", projectDirPath=projectDirectoryPath, rule="DGE")
  output:
    expand("{projectDirPath}csv/{{condition}}/{{condition}}.finishedDGE.txt", projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/15_dge.{{condition}}.txt", 3)
  script:
    "scripts/DGE.R"


rule createShinyApp:
  input:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    multiSampled = config["multiSampled"],
    cores = pOpt.numCores["shinyApp"],
    time=res.approxWalltime("shiny", sampleType, num_cells, additionalTime=pOpt.addTime["shinyApp"]),
    mem=res.approxRAM("shiny", sampleType, num_cells, additionalRAM=pOpt.addRAM["shinyApp"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="createShinyApp"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="createShinyApp")
  conda:
    pOpt.shinyEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="createShinyApp")
  output:
    expand("{projectDirPath}shinyApp/server.R", project=config["projectName"], projectDirPath=projectDirectoryPath),
    expand("{projectDirPath}shinyApp/ui.R", project=config["projectName"], projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/16_shiny.txt", 3)
  script:
    "scripts/createShinyApp.R"


rule multimodalFeaturePlotting:
  input:
    expand("{projectDirPath}outputs/{project}.markerDisc.rds", project=config["projectName"], projectDirPath=projectDirectoryPath)
  params:
    projectDirPath = projectDirectoryPath,
    project = config["projectName"],
    assayname = assayName,
    cores = pOpt.numCores["multimodalPlot"],
    time=res.approxWalltime("mmodalplot", sampleType, num_cells, additionalTime=pOpt.addTime["multimodalPlot"]),
    mem=res.approxRAM("mmodalplot", sampleType, num_cells, additionalRAM=pOpt.addRAM["multimodalPlot"]),
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="multimodalFeaturePlotting"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="multimodalFeaturePlotting")
  conda:
    pOpt.usualEnv
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="multimodalFeaturePlotting")
  output:
    expand("{projectDirPath}plots/finishedFeatures.txt", project=config["projectName"], projectDirPath=projectDirectoryPath)
#  benchmark:
#    repeat("benchmarks/" + config["projectName"] + "/17_mmPlot.txt", 3)
  script:
    "scripts/multimodalFeaturePlotting.R"


rule copyShinyAppInstructions:
  params:
    cores=1,
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="copyInstructions"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="copyInstructions")
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="copyInstructions")
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
    cores=1,
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingMTCutoff"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingMTCutoff")
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="missingMTCutoff")
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
    cores=1,
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingDoubletInfos"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingDoubletInfos")
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="missingDoubletInfos")
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
    cores=1,
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingConditions"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingConditions")
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="missingConditions")
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
    cores=1,
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingIntPCs"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingIntPCs")
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="missingIntPCs")
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
    cores=1,
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingResolutionsNeighPC"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingResolutionsNeighPC")
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="missingResolutionsNeighPC")
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
    cores=1,
    mem="50MB",
    time="0:01:00",
    error=expand("{projectDirPath}clusterLogs/{rule}.errors", projectDirPath=projectDirectoryPath, rule="missingChosenResolution"),
    output=expand("{projectDirPath}clusterLogs/{rule}.output", projectDirPath=projectDirectoryPath, rule="missingChosenResolution")
  log:
    expand("{projectDirPath}logs/{rule}.log", projectDirPath=projectDirectoryPath, rule="missingChosenResolution")
  output:
    temp(expand("{projectDirPath}workDirectory/chosenResolutionMissing.txt", projectDirPath=projectDirectoryPath))
  run:
    hf.raiseInputMissingException("errorMessage/chosenResolutionMissing.txt")
    shell("touch {output}")
