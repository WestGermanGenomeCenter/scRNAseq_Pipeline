#!/bin/bash

#load Snakemake and Conda for cluster execution on the HPC
module load Miniconda/3_snakemake
module load Snakemake/5.10.0

{exec_job}
