#!/bin/bash
yes | conda remove --name addmeta-env --all
yes | conda remove --name counting-env --all
yes | conda remove --name devtools --all
yes | conda remove --name doublet-env --all
yes | conda remove --name usual-env --all
yes | conda remove --name marker-env --all
yes | conda remove --name shiny-env --all

conda env create -f envs/addmeta.yml
conda env create -f envs/counting.yml
conda env create -f envs/devtools.yml
conda env create -f envs/doublet.yml
conda env create -f envs/env.yml
conda env create -f envs/marker.yml
conda env create -f envs/shiny.yml

conda activate addmeta-env
conda env export > envs/addmeta-spec.yml
conda deactivate
conda activate counting-env
conda env export > envs/counting-spec.yml
conda deactivate
conda activate devtools
conda env export > envs/devtools-spec.yml
conda deactivate
conda activate doublet-env
conda env export > envs/doublet-spec.yml
conda deactivate
conda activate usual-env
conda env export > envs/env-spec.yml
conda deactivate
conda activate marker-env
conda env export > envs/marker-spec.yml
conda deactivate
conda activate shiny-env
conda env export > envs/shiny-spec.yml
conda deactivate

python envs/env-cleaner.py