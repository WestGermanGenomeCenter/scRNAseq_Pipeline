import os
import hashlib


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
  folders = ["clusterLogs/", "csv/", "outputs/", "plots/", "shinyApp/", "workDirectory"]
  os.makedirs(projectDirectoryPath, exist_ok=True) #create directory if not exists
  for i in range(len(folders)):
      os.makedirs(projectDirectoryPath + folders[i], exist_ok=True) #create folder if not exists


def transform_sampleInputs(sampleInputs):
  samInputs = {}
  for i in sampleInputs[0].keys():
    samInputs[i] = []
  for i in sampleInputs:
    for j in i.keys():
      samInputs[j].append(i[j])
  return samInputs

  
def checkMinimalInputs(config):
  if not config["workDirectory"]:
    return False
  if not config["rawData"]:
    return False
  if not config["projectName"]:
    return False
  if not config["projectDirectoryPath"]:
    return False
  if not config["maxRAM"]:
    return False
  if config["multiSampled"] is None:
    return False
  if not config["numberOfCells"]:
    return False
  for i in config["sampleInputs"]:
    if i["name"] == None:
      return False
  return True
