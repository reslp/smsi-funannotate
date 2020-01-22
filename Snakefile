#singularity: "docker://reslp/funannotate:1.7.2"

configfile: "data/config.yaml"

import pandas as pd

sample_data = pd.read_table(config["samples"]).set_index("sample", drop=False)
#print(samples.assembly_path.tolist())

def get_assembly_path(wildcards):
# this is to get the assembly path information for the sample from the CSV file
	return sample_data.loc[wildcards.sample, ["assembly_path"]].to_list()


rule all:
	input:
		expand("results/output_{name}.fas", name=sample_data.index.tolist())

rule clean:
	input:
		assembly = get_assembly_path
	output:
		"results/output_{sample}.fas"
	shell:
		"""
		#touch {output}
		cp {input.assembly} {output}
		"""