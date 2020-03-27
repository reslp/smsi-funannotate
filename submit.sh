#!/bin/bash

module load singularity/3.5.2-gcc-9.1.0-fp2564h 
echo "the submission script is in testing mode! Check the snakemake command inside before moving into production!"


snakemake --use-singularity --singularity-args "-B $(pwd)/data/eggnogdb:/data/eggnogdb -B $(pwd)/data/database:/data/database -B $(pwd)/data/external:/data/external -B $(pwd)/data/RepeatMaskerLibraries:/software/RepeatMasker/Libraries" --jobs 1001 --cluster-config data/cluster_config-vsc4.yaml --cluster '$(pwd)/bin/immediate_submit.py {dependencies}' --immediate-submit -pr --notemp --latency-wait 600 $1

