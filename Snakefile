singularity: "docker://reslp/funannotate:1.7.2"

configfile: "data/config.yaml"

import pandas as pd
import os
WD=os.getcwd()

sample_data = pd.read_table(config["samples"]).set_index("sample", drop=False)
#print(samples.assembly_path.tolist())

def get_assembly_path(wildcards):
# this is to get the assembly path information for the sample from the CSV file
	return sample_data.loc[wildcards.sample, ["assembly_path"]].to_list()

def get_contig_prefix(wildcards):
	return sample_data.loc[wildcards.sample, ["contig_prefix"]].to_list()

rule all:
	input:
		expand("results/{name}/{name}_sorted.fas", name=sample_data.index.tolist())

rule clean:
	input:
		assembly = get_assembly_path
	output:
		"results/{sample}/{sample}_cleaned.fas"
	log:
		"log/{sample}_clean.log"
	params:
		folder = "{sample}"
	shell:
		"""
		if [[ ! -d results/{params.folder} ]]
		then
			mkdir results/{params.folder}
		fi
		cd results/{params.folder}
		funannotate clean -i ../../{input.assembly} -o ../../{output} 2> ../../{log}
		"""

rule sort:
	input:
		assembly = rules.clean.output
	output:
		"results/{samples}/{sample}_sorted.fas"
	params:
		folder = "{sample}",
		contig_prefix = get_contig_prefix
	shell:
		"""
                cd results/{params.folder}
                funannotate sort -i ../../{input.assembly} -o ../../{output} -b {params.contig_prefix}
		"""

