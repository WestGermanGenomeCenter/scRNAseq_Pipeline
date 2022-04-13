import time
import datetime
import math

memory = {
    "nm": {
        "meta": 0.00025,
        "mt1": 0.00015,
        "mt2": 0.00015,
        "dbElb": 0.00045,
        "dbRem": 0.00075,
        "merge": 0.00055,
        "sct": 0.00087,
        "integration": 0.00180,
        "UMAP": 0.00030,
        "tCluster": 0.00028,
        "cCluster": 0.00030,
        "mmodalAssay": 0.0,
        "marker": 0.00089,
        "count": 0.00115,
        "DGE": 0.00075,
        "shiny": 0.00099,
        "mmodalplot": 0.0,
        "demultiplex": 0
    },
    "ns": {
        "meta": 0.00046,
        "mt1": 0.00026, 
        "mt2": 0.00026, 
        "dbElb": 0.00139, 
        "dbRem": 0.00256, 
        "merge": 0.00046, 
        "sct": 0.00163, 
        "integration": 0.00046, 
        "UMAP": 0.00046, 
        "tCluster": 0.00046,
        "cCluster": 0.00046,
        "mmodalAssay": 0.0, 
        "marker": 0.00139, 
        "count": 0.00070, 
        "DGE": 0, 
        "shiny": 0.00048,
        "mmodalplot": 0.0,
        "demultiplex": 0
    },
    "mm": {
        "meta": 0.00026,
        "mt1": 0.00011, 
        "mt2": 0.00017, 
        "dbElb": 0.00039, 
        "dbRem": 0.00059, 
        "merge": 0.00055, 
        "sct": 0.00077, 
        "integration": 0.00150, 
        "UMAP": 0.00029, 
        "tCluster": 0.00027,
        "cCluster": 0.00026,
        "mmodalAssay": 0.00029, 
        "marker": 0.00089, 
        "count": 0.00088, 
        "DGE": 0.00062, 
        "shiny": 0.00060,
        "mmodalplot": 0.00045,
        "demultiplex": 0.0
    },
    "ms": {
        "meta": 0.00040,
        "mt1": 0.00040, 
        "mt2": 0.00078, 
        "dbElb": 0.00116, 
        "dbRem": 0.00192, 
        "merge": 0.00040, 
        "sct": 0.00154, 
        "integration": 0.00078, 
        "UMAP": 0.00078, 
        "tCluster": 0.00078,
        "cCluster": 0.00078,
        "mmodalAssay": 0.00078, 
        "marker": 0.00116, 
        "count": 0.00078, 
        "DGE": 0, 
        "shiny": 0.00078,
        "mmodalplot": 0.00078,
        "demultiplex": 0.0
    },
    "HTO": {
        "meta": 0,
        "mt1": 0.00076, 
        "mt2": 0.00076, 
        "dbElb": 0.00076, 
        "dbRem": 0, 
        "merge": 0.00114, 
        "sct": 0.00152, 
        "integration": 0.00267, 
        "UMAP": 0.00076, 
        "tCluster": 0.00076,
        "cCluster": 0.00114,
        "mmodalAssay": 0, 
        "marker": 0.00190, 
        "count": 0.00152, 
        "DGE": 0, 
        "shiny": 0.00100,
        "mmodalplot": 0.00114,
        "demultiplex": 0.00114
    }
}

w_time = {
    "nm": {
        "meta": 0.01169,
        "mt1": 0.00794, 
        "mt2": 0.00794, 
        "dbElb": 0.04086, 
        "dbRem": 0.30000, 
        "merge": 0.01605, 
        "sct": 0.03817, 
        "integration": 0.20356, 
        "UMAP": 0.02224, 
        "tCluster": 0.00607,
        "cCluster": 0.01967, 
        "mmodalAssay": 0.0, 
        "marker": 0.25456, 
        "count": 0.00794, 
        "DGE": 0.03258, 
        "shiny": 0.00741,
        "mmodalplot": 0.0,
        "demultiplex": 0
    },
    "ns": {
        "meta": 0.01743,
        "mt1": 0.01325, 
        "mt2": 0.01302, 
        "dbElb": 0.07183, 
        "dbRem": 0.45397, 
        "merge": 0.02441, 
        "sct": 0.07020, 
        "integration": 0.02510, 
        "UMAP": 0.03045, 
        "tCluster": 0.01046,
        "cCluster": 0.02487, 
        "mmodalAssay": 0.0,
        "marker": 0.12064, 
        "count": 0.03184, 
        "DGE": 0, 
        "shiny": 0.02092,
        "mmodalplot": 0.0,
        "demultiplex": 0
    },
    "mm": {
        "meta": 0.02140,
        "mt1": 0.00849, 
        "mt2": 0.00902, 
        "dbElb": 0.04979, 
        "dbRem": 0.23961, 
        "merge": 0.01793, 
        "sct": 0.03996, 
        "integration": 0.18380, 
        "UMAP": 0.02158, 
        "tCluster": 0.00366,
        "cCluster": 0.01952,
        "mmodalAssay": 0.01849, 
        "marker": 0.24493, 
        "count": 0.00500, 
        "DGE": 0.00786, 
        "shiny": 0.00553,
        "mmodalplot": 0.00326,
        "demultiplex": 0.03429
    },
    "ms": {
        "meta": 0.02629,
        "mt1": 0.01714, 
        "mt2": 0.01905, 
        "dbElb": 0.07086, 
        "dbRem": 0.40876, 
        "merge": 0.02629, 
        "sct": 0.06514, 
        "integration": 0.02857, 
        "UMAP": 0.03200, 
        "tCluster": 0.02476,
        "cCluster": 0.02743,
        "mmodalAssay": 0.02629, 
        "marker": 0.16927, 
        "count": 0.01143, 
        "DGE": 0, 
        "shiny": 0.01600,
        "mmodalplot": 0.00076,
        "demultiplex": 0.0
    },
    "HTO": {
        "meta": 0,
        "mt1": 0.03162, 
        "mt2": 0.03200, 
        "dbElb": 0.07771, 
        "dbRem": 0, 
        "merge": 0.03162, 
        "sct": 0.07733, 
        "integration": 0.11429, 
        "UMAP": 0.03771, 
        "tCluster": 0.03162,
        "cCluster": 0.06857,
        "mmodalAssay": 0, 
        "marker": 0.14286, 
        "count": 0.01714, 
        "DGE": 0, 
        "shiny": 0.02286,
        "mmodalplot": 0.01714,
        "demultiplex": 0.08000
    }
}

def approxRAM(rule, typ, num_cells, additionalRAM=0):
    #print(math.ceil(memory[typ][rule]*num_cells))
    approxMem = math.ceil(memory[typ][rule]*num_cells)+additionalRAM
    if approxMem == 0:
        approxMem = 1
    return str(approxMem) + "GB"

def approxWalltime(rule, typ, num_cells, additionalTime=0):
    #print(rule + " used: "+ time.strftime("%H:%M:%S", time.gmtime(math.ceil(w_time[typ][rule]*num_cells))), end=", ")
    apprTime = math.ceil(w_time[typ][rule]*num_cells)+additionalTime
    if apprTime < 120:
        apprTime = 300
    #print(apprTime)
    #print(str(datetime.timedelta(seconds=apprTime)))
    return str(datetime.timedelta(seconds=apprTime))
