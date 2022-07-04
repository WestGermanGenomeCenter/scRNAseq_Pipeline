
files = ["addmeta-spec.yml", "counting-spec.yml", "devtools-spec.yml", "doublet-spec.yml", "env-spec.yml", "marker-spec.yml", "shiny-spec.yml"]

bioconductorSeuratHelper = "\
  - bioconductor-biobase=2.54.0\n\
  - bioconductor-biocgenerics=0.40.0\n\
  - bioconductor-delayedarray=0.20.0\n\
  - bioconductor-iranges=2.28.0\n\
  - bioconductor-matrixgenerics=1.6.0\n\
  - bioconductor-s4vectors=0.32.0\n\
  - bioconductor-xvector=0.34.0\n\
  - bioconductor-zlibbioc=1.40.0\n\
"


for filename in files:
    textlines = []
    with open("envs/netAccess/" + filename) as f:
        textlines = f.readlines()
    index = 0
    while("dependencies:" not in textlines[index]):
        if "- conda-forge" in textlines[index]:
            textlines[index] = ""
        elif "- bioconda" in textlines[index]:
            textlines[index] = ""
        if "defaults" in textlines[index]:
            textlines[index] = "  - nodefaults\n"
        index += 1
    for i in range(index, len(textlines)):
        j = textlines[i].rfind("=")
        textlines[i] = textlines[i][0:j] + "\n"
    if "prefix" in textlines[-1]:
        textlines[-1] = ""
    with open("envs/HHU_HPC/" + filename, "w") as f:
        f.write("".join(textlines))