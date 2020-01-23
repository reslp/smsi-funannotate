#!/usr/bin/env python3

import os
import sys

from snakemake.utils import read_job_properties

# last command-line argument is the job script
jobscript = sys.argv[-1]

# all other command-line arguments are the dependencies
dependencies = set(sys.argv[1:-1])

# parse the job script for the job properties that are encoded by snakemake within
# this also includes the information contained in the cluster-config file as job_properties["cluster"]
job_properties = read_job_properties(jobscript)

# create list with command line arguments
cmdline = ["sbatch"]

# create string for slurm submit options for rule
slurm_args = "--partition={partition} --qos={qos} --mem={mem} --ntasks={ntasks} --ntasks-per-node={ntasks-per-node} --time={time} --hint={hint} --output={output} --output={output} -N {N} -J {J}".format(**job_properties["cluster"])
cmdline.append(slurm_args)

# now work on dependencies
if dependencies:
    cmdline.append("--dependency")
    # only keep numbers (which are the jobids) in dependencies list. this is necessary because slurm returns more than the jobid. For other schedulers this could be different!
    dependencies = [x for x in dependencies if x.isdigit()]
    cmdline.append("afterok:" + ",".join(dependencies))

cmdline.append(jobscript)

#now write final commandback to the system
os.system(" ".join(cmdline))

