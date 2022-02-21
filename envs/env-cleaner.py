
files = ["counting-spec.yml", "devtools-spec.yml", "doublet-spec.yml", "env-spec.yml", "marker-spec.yml", "shiny-spec.yml"]

for filename in files:
    textlines = []
    with open(filename) as f:
        textlines = f.readlines()
    index = 0
    while("dependencies:" not in textlines[index]):
        index += 1
    for i in range(index, len(textlines)):
        j = textlines[i].rfind("=")
        textlines = textlines[0:j]
    with open(filename, "w") as f:
        f.write("\n".join(textlines))