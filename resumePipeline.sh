#!/bin/bash

if [[ $# -eq 0 ]]
then
        echo "Please use the configfile you want as argument like 'bash clusterExecution.sh path/to/config.yaml'"
else
        echo "$#"
        module load Miniconda/3_snakemake
        module load Snakemake/5.10.0

        snakemake --unlock --configfile $1 --use-conda --cluster-config cluster/cluster.json --cluster "qsub -A {cluster.account} -q {cluster.queue} -l select={cluster.nodes}:ncpus={cluster.ppn}:mem={params.mem} -l walltime={params.time} -o {params.output} -e {params.error}" --jobs 20 -p --keep-going --jobscript cluster/jobscript.sh --latency-wait 180
fi
