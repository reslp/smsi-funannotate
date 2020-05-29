singularity: "docker://reslp/funannotate:1.7.2"

configfile: "data/config.yaml"

import pandas as pd
import os

# useful variable definition:
WD=os.getcwd()
email="philipp.resl@uni-graz.at"

sample_data = pd.read_table(config["samples"]).set_index("sample", drop=False)
def get_assembly_path(wildcards):
# this is to get the assembly path information for the sample from the CSV file
	return sample_data.loc[wildcards.sample, ["assembly_path"]].to_list()

def get_contig_prefix(wildcards):
	return sample_data.loc[wildcards.sample, ["contig_prefix"]].to_list()

def get_premasked_state(wildcards):
	return sample_data.loc[wildcards.sample, ["premasked"]].to_list()[-1]

def get_all_samples(wildcards):
	sam = sample_data["contig_prefix"].tolist()
	sam = [sample+"s" for sample in sam].join()
	return sam 

rule all:
	input:
		#expand("results/{name}/{name}_remote.done", name=sample_data.index.tolist()),
		#expand("results/{name}/{name}_iprscan.done", name=sample_data.index.tolist()),
		#expand("results/{name}/{name}_eggnog.done", name=sample_data.index.tolist()),
		expand("results/{name}/{name}_tarpredict.done", name=sample_data.index.tolist()),
		expand("results/{name}/{name}_predict.done", name=sample_data.index.tolist()),
		expand("results/{name}/{name}_annotate.done", name=sample_data.index.tolist()),
		"results/funannotate_compare.done"

include: "rules/funannotate.smk"
