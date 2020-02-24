#!/bin/bash

module load singularity/3.5.2-gcc-9.1.0-fp2564h 

snakemake --use-singularity --singularity-args "-B $(pwd)/data/eggnogdb:/data/eggnogdb -B $(pwd)/data/database:/data/database -B $(pwd)/data/external:/data/external" --jobs 1001 --cluster-config data/cluster_config-vsc4.yaml --cluster '$(pwd)/bin/immediate_submit.py {dependencies}' --immediate-submit -pr --notemp --latency-wait 600 $1

