import matplotlib.pyplot as plt
import pandas as pd

data = pd.DataFrame([["165_1core", "integration", 6200], ["215_1core", "integration", 640],
                     ["165_2cores", "integration", 3800], ["215_2cores", "integration", 510],
                     ["165_1core", "markerDisc", 10700], ["215_1core", "markerDisc", 850],
                     ["165_2cores", "markerDisc", 6900], ["215_2cores", "markerDisc", 610],
                     ["165_1core", "DGE", 460], ["215_1core", "DGE", 160],
                     ["165_2cores", "DGE", -1], ["215_2cores", "DGE", 150]],
                    columns=["sample_numCores", "rule", "time_in_s"]
                    )
