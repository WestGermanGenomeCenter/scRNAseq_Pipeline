import time
import datetime
import math
import yaml
import os
import snakemakeScripts.cluster.plotRessources as pltRes

if not (os.path.exists("memoryGradient.yaml") and os.path.exists("walltimeGradient.yaml")):
    pltRes.createResourceYAML(False)

with open("memoryGradient.yaml", "r") as RAM:
    memory = yaml.safe_load(RAM)

with open("walltimeGradient.yaml", "r") as TIME:
    w_time = yaml.safe_load(TIME)

def approxRAM(rule, typ, num_cells, additionalRAM=0):
    #print(math.ceil(memory[typ][rule]*num_cells))
    approxMem = math.ceil(memory[typ][rule]*num_cells)+additionalRAM
    if approxMem == 0:
        approxMem = 1
    return str(approxMem) + "GB"

def approxWalltime(rule, typ, num_cells, additionalTime=0):
    #print(rule + " used: "+ time.strftime("%H:%M:%S", time.gmtime(math.ceil(w_time[typ][rule]*num_cells))), end=", ")
    apprTime = math.ceil(w_time[typ][rule]*num_cells)+additionalTime
    if apprTime < 300:
        apprTime = 300
    #print(apprTime)
    #print(str(datetime.timedelta(seconds=apprTime)))
    return str(datetime.timedelta(seconds=apprTime))
