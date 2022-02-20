conda env create -f envs/counting.yml
conda env create -f envs/devtools.yml
conda env create -f envs/doublet.yml
conda env create -f envs/env.yml
conda env create -f envs/marker.yml
conda env create -f envs/shiny.yml

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