import time
import math
import os
import matplotlib.pyplot as plt


def createNsavePlot(rule, func, typ, xvalues, yvalues, ylabel, leg, rulenumber):
    plt.figure()
    plot, = plt.plot(xvalues, yvalues, "o", label=typ)
    #plotfolder = "/plots/ns/"
    #plotfolder = "/plots/ms/"
    #plotfolder = "/plots/nm/"
    #plotfolder = "/plots/mm/"
    plotfolder = "/plots/HTO/"
    plt.title(rule + ": " + func)
    plt.xlabel("Cells per samples (average)")
    plt.ylabel(ylabel)
    legend = plt.legend(handles=[plot], title=leg, bbox_to_anchor=(1.05, 1))
    plt.savefig(os.path.dirname(os.path.abspath(__file__)) + plotfolder + rulenumber + "_" + rule + "_" + typ +
                ".png", bbox_extra_artists=(legend,), bbox_inches="tight")
    plt.close()


rules = ["meta", "mt1", "mt2", "dbElbow", "dbr", "merge", "sct", "integration", "UMAP", "testCluster", "chosenCluster",
         "marker", "count", "DGE", "shiny", "mmAssay", "mmPlot", "demux"]

times = {
    "165": {
        "RAM": [12, 4.8, 7.3, 20, 30, 27, 35, 72, 14, 13, 13, 45, 40, 29, 24, 0, 0, 0],
        "walltime": [397, 272, 285, 1366, 6932, 572, 1211, 8962, 979, 245, 866, 9148, 184, 230, 182, 0, 0, 0],
        "multi-sample": True,
        "multi-modal": False,
        "HTO": False,
        "cells": 66040
    },
    "215": {
        "RAM": [1.66, 0.9, 1.4, 3, 5, 5, 6.3, 17.6, 2.5, 1.6, 2.5, 7.6, 2.5, 5.5, 3.1, 0, 0, 0],
        "walltime": [60, 42, 42, 300, 2038, 100, 285, 460, 135, 45, 119, 1347, 70, 287, 65, 0, 0, 0],
        "multi-sample": True,
        "multi-modal": False,
        "HTO": False,
        "cells": 13228
    },
    "143_165": {
        "RAM": [1, 0.78, 1.5, 5, 8.6, 1.1, 4.6, 1.7, 1.2, 1.3, 1.2, 4.5, 2.1, 0, 1.3, 0, 0, 0],
        "walltime": [60, 40, 41, 233, 1436, 67, 221, 70, 95, 42, 70, 557, 20, 0, 41, 0, 0, 0],
        "multi-sample": False,
        "multi-modal": False,
        "HTO": False,
        "cells": 7757
    },
    "143_215": {
        "RAM": [1, 0.79, 1.4, 3.4, 8.2, 1, 4, 1.6, 1, 1.1, 1, 4, 2, 0, 1.2, 0, 0, 0],
        "walltime": [57, 40, 40, 212, 1290, 70, 211, 67, 88, 22, 68, 571, 30, 0, 40, 0, 0, 0],
        "multi-sample": False,
        "multi-modal": False,
        "HTO": False,
        "cells": 7757
    },
    "244_165": {
        "RAM": [10, 5, 5, 23, 35, 26, 42, 78, 14, 13, 13, 41, 51, 35, 30, 0, 0, 0],
        "walltime": [520, 353, 353, 1817, 10394, 714, 1698, 5918, 941, 180, 833, 11323, 205, 380, 300, 0, 0, 0],
        "multi-sample": True,
        "multi-modal": False,
        "HTO": False,
        "cells": 66723
    },
    "244_215": {
        "RAM": [13, 5, 5, 23, 35, 26, 44, 78, 15, 14, 14, 45, 55, 33, 32, 0, 0, 0],
        "walltime": [372, 263, 264, 1293, 6836, 582, 1272, 4408, 776, 270, 694, 10668, 310, 245, 232, 0, 0, 0],
        "multi-sample": True,
        "multi-modal": False,
        "HTO": False,
        "cells": 66723
    },
    "305_165": {
        "RAM": [1.2, 0.73, 0.72, 4.5, 7.8, 1.1, 5.3, 1.38, 1.2, 1.1, 1.2, 3.76, 2, 0, 1.3, 0, 0, 0],
        "walltime": [49, 37, 37, 205, 1291, 70, 190, 67, 82, 30, 70, 298, 16, 0, 40, 0, 0, 0],
        "multi-sample": False,
        "multi-modal": False,
        "HTO": False,
        "cells": 4302
    },
    "305_215": {
        "RAM": [1.2, 0.73, 0.73, 4.6, 7.8, 1.2, 4.8, 1.5, 1.3, 1.2, 1.2, 3.9, 2.1, 0, 1.5, 0, 0, 0],
        "walltime": [50, 38, 37, 206, 1302, 70, 201, 72, 87, 22, 71, 346, 120, 16, 60, 0, 0, 0],
        "multi-sample": False,
        "multi-modal": False,
        "HTO": False,
        "cells": 4302
    },
    "165_mm": {
        "RAM": [13, 5, 8, 20, 30, 27, 39, 76, 14, 13.2, 13, 45, 44, 31, 26, 14, 23, 0],
        "walltime": [942, 374, 397, 2192, 10549, 789, 1759, 8092, 950, 161, 859, 10783, 220, 346, 243, 814, 143, 0],
        "multi-sample": True,
        "multi-modal": True,
        "HTO": False,
        "cells": 66040
    },
    "143_165_mm": {
        "RAM": [1.5, 0.8, 1.3, 5.2, 9, 1.2, 5, 1.7, 1.3, 1.2, 1.1, 4.6, 2.1, 0, 1.4, 1.3, 1.1, 0],
        "walltime": [117, 49, 52, 310, 1916, 78, 277, 86, 113, 16, 93, 836, 16, 0, 45, 80, 24, 0],
        "multi-sample": False,
        "multi-modal": True,
        "HTO": False,
        "cells": 7757
    },
    "143_215_mm": {
        "RAM": [1.5, 0.8, 1.3, 3.5, 9, 1.1, 3.3, 1.7, 1.2, 1.2, 1.2, 4.4, 1.9, 0, 2.5, 1.3, 1.2, 0],
        "walltime": [102, 57, 47, 285, 1635, 66, 264, 97, 104, 15, 83, 875, 16, 0, 58, 88, 23, 0],
        "multi-sample": False,
        "multi-modal": True,
        "HTO": False,
        "cells": 7757
    },
    "284": {
        "RAM": [0, 0.9, 1.1, 1.5, 0, 2,2.2, 3, 1, 1, 1, 3.2, 3, 0, 1.5, 0, 1.5, 2.2],
        "walltime": [0, 55, 50, 130, 0, 55, 125, 177, 66, 55, 120, 230, 20, 0, 35, 0, 30, 133],
        "multi-sample": False,
        "multi-modal": True,
        "HTO": True,
        "cells": 2625
    },
    "284_165": {
        "RAM": [0, 0.9, 1, 1.5, 0, 2, 2.4, 5, 1, 1, 2, 3, 1.5, 0, 1.4, 0, 1.6, 2.2],
        "walltime": [0, 50, 55, 136, 0, 55, 135, 200, 62, 55, 120, 250, 30, 0, 40, 0, 25, 130],
        "multi-sample": False,
        "multi-modal": True,
        "HTO": True,
        "cells": 2625
    },
    "284_215": {
        "RAM": [0, 0.9, 1, 1.5, 0, 1.9, 2.3, 3.1, 1, 1, 2, 3.3, 1.5, 0, 1.5, 0, 1.5, 2.2],
        "walltime": [0, 50, 55, 135, 0, 55, 130, 180, 60, 55, 120, 250, 20, 0, 35, 0, 25, 140],
        "multi-sample": False,
        "multi-modal": True,
        "HTO": True,
        "cells": 2625
    }
}

for i in range(len(rules)):
    #print(rules[i])
    xaxis = []
    walltime = []
    ram = []
    rlegend = ""
    wlegend = ""
    max_wGradient = 0
    max_rGradient = 0
    for j in times.keys():
        #if not times[j]["multi-sample"] and not times[j]["multi-modal"]: #ns
        #if not times[j]["multi-sample"] and times[j]["multi-modal"]: #ms
        #if times[j]["multi-sample"] and not times[j]["multi-modal"]: #nm
        #if times[j]["multi-sample"] and times[j]["multi-modal"]: #mm
        if times[j]["HTO"]:
            rlegend += j + " -> "
            wlegend += j + " -> "
            wlegend += "Walltime: " + time.strftime("%H:%M:%S", time.gmtime(math.ceil(1.5 * times[j]["walltime"][i]))) +"\n"
            rlegend += "RAM: " + str(math.ceil(1.3 * times[j]["RAM"][i])) + "\n"
            newWalltime = math.ceil(1.5 * times[j]["walltime"][i])
            newRam = math.ceil(1.3 * times[j]["RAM"][i])
            walltime.append(newWalltime)
            ram.append(newRam)
            xaxis.append(times[j]["cells"])
            if newWalltime/float(times[j]["cells"]) > max_wGradient:
                #print(j)
                max_wGradient = newWalltime/float(times[j]["cells"])
            if newRam/float(times[j]["cells"]) > max_rGradient:
                #print(j)
                max_rGradient = newRam/float(times[j]["cells"])
    createNsavePlot(rules[i], ("f(x)="+"{:.5f}".format(max_rGradient)+"x"), "ram", xaxis, ram, "RAM in GB", rlegend, str(i))
    createNsavePlot(rules[i], ("f(x)="+"{:.5f}".format(max_wGradient)+"x"), "walltime", xaxis, walltime, "RAM in GB", wlegend, str(i))
    #print()