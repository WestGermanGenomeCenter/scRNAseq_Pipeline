import os
import hashlib
import itertools
from collections import OrderedDict



class MissingInputError(Exception):
    """Exception raised for missing input in configfile

    Attributes:
        message -- which input is mising
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def raiseInputMissingException(filename):
    msg = None
    with open(filename) as f:
        msg = f.readlines()
    raise MissingInputError(message=msg)


# taken from the conda.py file from snakemake to replicating the hashes for the environment folders
def findHash(env_yml):
  dir = os.path.realpath(".snakemake/conda")
  content = None
  with open(env_yml, "rb") as f:
    content = f.read()
  md5hash = hashlib.md5()
  md5hash.update(dir.encode())
  md5hash.update(content)
  return ".snakemake/conda/" + md5hash.hexdigest()[:8]


def createDirectoriesIfNotExists(projectDirectoryPath):
  folders = ["clusterLogs/", "csv/", "logs/", "outputs/", "plots/", "shinyApp/", "workDirectory"]
  os.makedirs(projectDirectoryPath, exist_ok=True) #create directory if not exists
  for i in range(len(folders)):
      os.makedirs(projectDirectoryPath + folders[i], exist_ok=True) #create folder if not exists


def transform_sampleInputs(sampleInputs):
  samples = OrderedDict()
  otherMetaData = []
  for i in sampleInputs:
    samples[i['name']] = i
  for i in sampleInputs:
    otherMetaData.append(i["otherMetaData"])
  return samples, otherMetaData


def determineSampleTypeAndAssayName(config):
  sampleType = "ns"
  assayName = "ADT"
  if config["multiSampled"] and config["multimodal"]:
    sampleType = "mm"
  elif config["multiSampled"] and not config["multimodal"]:
    sampleType = "nm"
  elif not config["multiSampled"] and config["multimodal"]:
    sampleType = "ms"
  if(config["HTO"]):
    assayName = "HTO"
    sampleType = "HTO"
  return sampleType, assayName

def checkMinimalInputs(config):
  #general parameters
  if config["HHU_HPC"] is None:
    return False
  if not config["maxRAM"]:
    return False
  #sample type
  if config["multiSampled"] is None:
    return False
  if config["multimodal"] is None:
    return False
  if config["HTO"] is None:
    return False
  #directories and inputs
  if not config["workDirectory"]:
    return False
  if not config["rawData"]:
    return False
  if not config["projectDirectoryPath"]:
    return False
  #sample information
  if not config["projectName"]:
    return False
  if not config["numberOfCells"]:
    return False
  if not config["mtPattern"]:
    return False
  if not config["otherMetaName"] or not type(config["otherMetaName"]) is list:
    return False
  for i in config["sampleInputs"]:
    if i["name"] == None:
      return False
    if not i["otherMetaData"]:
      return False
    if not i["expectedPercentDoublets"]:
      return False
  return True


def createMultiSampleInput(path, folder, samples, ending, project=None):
  ioputs = []
  for i in samples:
    ioputs.append(path + folder + str(i) + ending)
  if project: #Project is only not none at the very beginning of an not HTO, multimodal project
    ioputs.append(path + folder + project + ".rawData.rds")
  return ioputs

def createCombinations(conditions, combiOnly=True):
  output = []
  for i in range(1, len(conditions)+1):
    output += itertools.combinations(conditions, i)
  for i in range(0, len(output)):
    output[i] = '_'.join(output[i])
  if(combiOnly):
    output = [o for o in output if o not in conditions]
  #print(output)
  return output

def createMultiMetaCountInput(path, folder, conditions, ending):
  ioputs = createCombinations(conditions, combiOnly=False)
  for i in range(0, len(ioputs)):
    ioputs[i] = path + folder + ioputs[i] + ending
  #print(ioputs)
  return ioputs
