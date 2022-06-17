import math
import os
import yaml
import matplotlib.pyplot as plt


rules = ["meta", "demultiplex", "mt1", "mt2", "dbElb", "dbRem", "merge",
        "sct", "integration", "UMAP", "tCluster", "cCluster", "mmodalAssay",
        "marker", "count", "DGE", "shiny", "mmodalplot"]

times = {
    "ns": {
        "143_165ns": {
            "RAM":      [1.1, -1, 0.85, 1.6, 5.1,  9.2, 1.1, 5.8, 1.5, 1.4, 1.3, 1.4, -1, 5.3, 1.5, -1, 1.5, -1],
            "walltime": [ 75, -1,   60,  60, 240, 1500,  80, 230,  75, 110,  40,  80, -1, 550,  30, -1,  50, -1],
            "cells": 7757
        },
        "143_215ns": {
            "RAM":      [1.1, -1, 0.9, 1.5, 3.5,  7.3, 1.2, 4.1, 1.5, 1.3, 1.3, 1.3, -1, 4.4, 1.1, -1, 0.9, -1],
            "walltime": [ 75, -1,  50,  50, 230, 1500,  75, 230,  80, 100,  30,  80, -1, 550,  40, -1,  60, -1],
            "cells": 7757
        },
        "305_165": {
            "RAM":      [1.4, -1, 0.8, 0.9,   6,    8, 1.2,   5, 1.5, 1.3, 1.3, 1.3, -1,   4, 1.5, -1, 2.4, -1],
            "walltime": [ 70, -1,  50,  50, 230, 1450,  90, 230,  90, 100,  40,  80, -1, 350,  50, -1,  70, -1],
            "cells": 4302
        },
        "305_215": {
            "RAM":      [1.4, -1, 0.8, 0.8,   6,    8, 1.3,   5, 1.6, 1.3, 1.4, 1.4, -1,   4, 1.5, -1, 1.4, -1],
            "walltime": [ 70, -1,  50,  50, 230, 1400, 100, 230,  90, 120,  30,  90, -1, 400,  40, -1,  50, -1],
            "cells": 4302
        }
    },
    "nm" : {
        "165nm": {
            "RAM":      [ 13, -1, 1.6, 2.7, 11.1,   19,  30,   40,   80,   14,  13,  14, -1,    45,  45,  32,  27, -1],
            "walltime": [800, -1, 110,  95,  410, 1900, 900, 1800, 6200, 1000, 230, 900, -1, 10800, 250, 500, 260, -1],
            "cells": 66040
        },
        "215": {
            "RAM":      [ 2, -1,  1,  1, 2.1,   3, 5.1, 6.7,  22, 2.7, 1.7, 2.8, -1,   8,  2, 6.5, 4.1, -1],
            "walltime": [90, -1, 30, 30, 115, 720, 160, 500, 700, 200,  50, 200, -1, 900, 60, 200,  60, -1],
            "cells": 13228
        },
        "244_165": {
            "RAM":      [ 11, -1,   2,   2,  16,   25,  27,   38,   80,   14,  13,  13, -1,    41,  46,  35,  27, -1],
            "walltime": [580, -1, 130, 130, 400, 2100, 800, 1780, 6500, 1000, 200, 900, -1, 14000, 300, 570, 300, -1],
            "cells": 66723
        },
        "244_215": {
            "RAM":      [ 11, -1,   2,   2,  16,   26,  27,   44,   83,   16,  14, 14.2, -1,    51,  47,  35,  29, -1],
            "walltime": [560, -1, 130, 130, 400, 2200, 950, 1900, 6900, 1100, 260, 1000, -1, 13800, 300, 600, 300, -1],
            "cells": 66723
        },
        "354": {
            "RAM":      [ 13, -1, 2.3, 2.1,  23,   35,  32,   40,   52,  13,  12,  12, -1,   46,  39,  40,  25, -1],
            "walltime": [600, -1, 150, 150, 600, 2450, 800, 1600, 2800, 880, 160, 800, -1, 4300, 250, 900, 250, -1],
            "cells": 46000
        },
        "354m": {
            "RAM":      [ 13, -1, 2.3, 4.1,  22,   33,  28,   47,   50,  13,  13,  13, -1,   47,  48,   44,  27, -1],
            "walltime": [550, -1, 150, 150, 600, 2450, 800, 1750, 2650, 900, 160, 850, -1, 5400, 300, 2350, 250, -1],
            "cells": 46000
        }
    },
    "ms": {
        "143_165ms": {
            "RAM":      [1.1, -1, 0.85, 1.6, 5.1,  9.2, 1.2, 5.8, 1.6, 1.4, 1.3, 1.4, 1.5, 4.7, 1.6, -1, 1.5, 1.5],
            "walltime": [ 75, -1,   60,  60, 240, 1500,  80, 240,  80, 110,  40,  80,  75, 550,  40, -1,  60,  60],
            "cells": 7757
        },
        "143_215ms": {
            "RAM":      [1.1, -1, 0.9, 1.5, 3.6,    8, 1.2, 4.1, 1.5, 1.3, 1.3, 1.3, 1.5, 4.2, 1.6, -1, 1.5, 1.5],
            "walltime": [ 80, -1,  60,  60, 230, 1500,  80, 230,  80, 100,  40,  80,  80, 550,  40, -1,  60,  60],
            "cells": 7757
        },
    },
    "mm": {
        "165mm": {
            "RAM":      [ 13, -1, 1.6, 3.2, 10.2,   18,  34,   38,   80,   15,  13,  14,  16,    46,  45,  33,  27,  23],
            "walltime": [750, -1,  90,  95,  400, 2200, 950, 1900, 6600, 1000, 230, 950, 900, 10500, 265, 550, 250, 220],
            "cells": 66040
        }
    }, 
    "HTO": {
        "116": {
            "RAM":      [-1,  18, 1.7, 1.7, 7.3, -1,  10,  16,   20, 5.5,   5,   5, -1,   21,  12,  17,  12, 10],
            "walltime": [-1, 600, 110, 100, 300, -1, 300, 800, 1400, 400, 100, 640, -1, 1750, 150, 900, 120, 95],
            "cells": 31000
        },
        "284": {
            "RAM":      [-1, 2.5, 0.6, 0.6, 0.6, -1, 1.8, 2.2, 5.5,  1,  1,  1, -1,   3, 1.5, 0, 1.5, 1.5],
            "walltime": [-1, 200,  30,  30,  50, -1,  60, 130, 200, 80, 30, 70, -1, 250,  40, 0,  50,  40],
            "cells": 2625
        },
        "284_165": {
            "RAM":      [-1,   2, 0.5, 0.6, 0.5, -1,  2, 2.2,   4,  1,  1,  1, -1,   3, 1.5, 2.5, 1.5, 1.5],
            "walltime": [-1, 150,  30,  30,  60, -1, 60, 140, 190, 70, 30, 60, -1, 230,  35,  80,  60,  30],
            "cells": 2625
        },
        "284_215": {
            "RAM":      [-1, 2.5, 0.5, 0.6, 0.5, -1,  2, 2.2,   6,  1,  1,  1, -1, 3.2, 1.5, 2.1, 1.5, 1.5],
            "walltime": [-1, 150,  40,  40,  40, -1, 60, 140, 200, 70, 30, 60, -1, 300,  40,  80,  60,  40],
            "cells": 2625
        }
    }
}


def createPlot(rule, sType, ruleNum, cells, RAM, TIME, memGrad, timeGrad):
    plotfolder = os.path.dirname(os.path.abspath(__file__)) + "/plots/" + sType + "/"
    os.makedirs(plotfolder, exist_ok=True)
    fig = plt.figure()

    timePlot = fig.add_subplot(2,1,1)
    t, = timePlot.plot(cells, TIME, "^", alpha=0.4)
    timePlot.set_title("Time=" + str(timeGrad) + "x")
    timePlot.set_xlabel("Cells per dataset")
    timePlot.set_ylabel("Walltime in seconds")

    ramPlot = fig.add_subplot(2,1,2)
    r, = ramPlot.plot(cells, RAM, "o", alpha=0.4)
    ramPlot.set_title("Memory=" + str(memGrad) + "y")
    ramPlot.set_xlabel("Cells per dataset")
    ramPlot.set_ylabel("Walltime in seconds")

    fig.tight_layout()
    legend = fig.legend(handles=[t, r], labels=["a", "b"], title="Hi", bbox_to_anchor=(1.05, 1))
    fig.savefig(plotfolder + str(ruleNum) + "_" + rule + ".png", bbox_extra_artists=(legend,)) 
    plt.close()


def createResourceYAML(createPlots=False):
    w_time = {"ns":{}, "nm":{}, "ms":{}, "mm":{}, "HTO":{}}
    memory = {"ns":{}, "nm":{}, "ms":{}, "mm":{}, "HTO":{}}
    for i in range(len(rules)):             #for every rule
        for j in times.keys():              #for every pipeline way
            gigaByte = []
            wSeconds = []
            numCells = []
            for k in times[j].keys():       #for every sample benchmarked
                gigaByte.append(times[j][k]["RAM"][i])
                wSeconds.append(times[j][k]["walltime"][i])
                numCells.append(times[j][k]["cells"])
            maxRAMGradient = 0
            maxSecGradient = 0
            for k in range(len(numCells)):
                w = math.ceil(1.5*float(wSeconds[k]))/float(numCells[k])
                if w > maxSecGradient:
                    maxSecGradient = w
                g = math.ceil(1.3*float(gigaByte[k]))/float(numCells[k])
                if g > maxRAMGradient:
                    maxRAMGradient = g
            maxRAMGradient = round(maxRAMGradient, 5)
            maxSecGradient = round(maxSecGradient, 5)
            w_time[j][rules[i]] = maxSecGradient
            memory[j][rules[i]] = maxRAMGradient
            if(createPlots):
                createPlot(rule=rules[i], ruleNum=i, sType=j, 
                        cells=numCells, RAM=gigaByte, TIME=wSeconds, 
                        memGrad=maxRAMGradient, timeGrad=maxSecGradient)


    with open("walltimeGradient.yaml", "w") as out:
        yaml.dump(w_time, out, default_flow_style=False)

    with open("memoryGradient.yaml", "w") as out:
        yaml.dump(memory, out, default_flow_style=False)

#createResourceYAML(False)