#!/bin/bash
yes | conda remove --name addmeta-env --all
yes | conda remove --name counting-env --all
yes | conda remove --name devtools --all
yes | conda remove --name doublet-env --all
yes | conda remove --name usual-env --all
yes | conda remove --name marker-env --all
yes | conda remove --name shiny-env --all

conda env create -f envs/netAccess/addmeta.yml
conda env create -f envs/netAccess/counting.yml
conda env create -f envs/netAccess/devtools.yml
conda env create -f envs/netAccess/doublet.yml
conda env create -f envs/netAccess/env.yml
conda env create -f envs/netAccess/marker.yml
conda env create -f envs/netAccess/shiny.yml

conda activate addmeta-env
conda env export > envs/HHU_HPC/addmeta-spec.yml
conda deactivate
conda activate counting-env
conda env export > envs/HHU_HPC/counting-spec.yml
conda deactivate
conda activate devtools
conda env export > envs/HHU_HPC/devtools-spec.yml
conda deactivate
conda activate doublet-env
conda env export > envs/HHU_HPC/doublet-spec.yml
conda deactivate
conda activate usual-env
conda env export > envs/HHU_HPC/env-spec.yml
conda deactivate
conda activate marker-env
conda env export > envs/HHU_HPC/marker-spec.yml
conda deactivate
conda activate shiny-env
conda env export > envs/HHU_HPC/shiny-spec.yml
conda deactivate

python envs/env-cleaner.py