singularity: "docker://reslp/funannotate:latest"

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

rule all:
	input:
		#expand("results/{name}/{name}_remote.done", name=sample_data.index.tolist()),
		#expand("results/{name}/{name}_iprscan.done", name=sample_data.index.tolist()),
		#expand("results/{name}/{name}_eggnog.done", name=sample_data.index.tolist()),
		expand("results/{name}/{name}_tarpredict.done", name=sample_data.index.tolist()),
		expand("results/{name}/{name}_annotate.done", name=sample_data.index.tolist())

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
		funannotate clean -i ../../{input.assembly} -o ../../{output} --minlen 5000 2> ../../{log}
		"""

rule sort:
	input:
		assembly = rules.clean.output
	output:
		"results/{sample}/{sample}_sorted.fas"
	params:
		folder = "{sample}",
		contig_prefix = get_contig_prefix
	shell:
		"""
		cd results/{params.folder}
		funannotate sort -i ../../{input.assembly} -o ../../{output} -b {params.contig_prefix}
		"""

rule mask:
	input:
		assembly = rules.sort.output
	output:
		"results/{sample}/{sample}_masked.fas"
	params:
		folder = "{sample}",
		premasked = get_premasked_state
	shell:
		"""
		if [[ "{params.premasked}" == "yes" ]]; then
			cp {input.assembly} {output}
		else
			cd results/{params.folder}
			funannotate mask -i ../../{input.assembly} -o ../../{output} -m repeatmasker
		fi
		"""	


rule predict:
	input:
		assembly = rules.mask.output 
	output:
		"results/{sample}/{sample}_predict.done"
	params:
		folder = "{sample}",
		pred_folder = get_contig_prefix,
		sample_name = "{sample}",
	log:
		"log/{sample}_predict.log"
	threads: 16
	shell:
		"""
		if [[ ! -d results/{params.folder} ]]
		then
			mkdir results/{params.folder}
		fi
		cd results/{params.folder}
		funannotate predict -i ../../{input.assembly} -o {params.pred_folder}_preds -s {params.sample_name} --name {params.pred_folder}_pred --optimize_augustus --cpus {threads} --busco_db metazoa --organism other --busco_seed_species schistosoma --ploidy 2 >& ../../{log}
		touch ../../{output}
		""" 

rule tarpredict:
	input:
		{rules.predict.output}
	output:
		"results/{sample}/{sample}_tarpredict.done"
	params:
		pred_folder = get_contig_prefix,
		folder = "{sample}"
	shell:
		"""
		cd results/{params.folder}/{params.pred_folder}_preds/predict_misc
		tar -cvf EVM_busco.tar EVM_busco && rm -r EVM_busco
		tar -cvf busco.tar busco && rm -r busco
		tar -cvf genemark.tar genemark && rm -r genemark
		tar -cvf busco_proteins.tar busco_proteins && rm -r busco_proteins
		tar -cvf EVM.tar EVM && rm -r EVM
		touch ../../../../{output}
		"""		
rule iprscan:
	input:
		rules.predict.output

	output:
		"results/{sample}/{sample}_iprscan.done"
	params:
		folder="{sample}",
		pred_folder=get_contig_prefix
	log:
		"log/{sample}_ipscan.log"
	shell:
		"""
		cd results/{params.folder}
		funannotate iprscan --iprscan_path /data/external/interproscan-5.39-77.0/interproscan.sh -i ../../results/{params.folder}/{params.pred_folder}_preds -m local -c 2 >& ../../{log}
		touch ../../{output}
		"""
rule remote:
	input:
		rules.predict.output
	output:
		"results/{sample}/{sample}_remote.done"
	params:
		folder="{sample}",
		pred_folder=get_contig_prefix
	log:
		"log/{sample}_remote.log"
	shell:
		"""
		cd results/{params.folder}
		funannotate remote -i {params.pred_folder}_preds -m phobius -e philipp.resl@uni-graz.at >& ../../{log}
		touch ../../{output}
		"""
rule eggnog:
	input:
		rules.predict.output
	output:
		"results/{sample}/{sample}_eggnog.done"
	params:
		folder="{sample}",
		pred_folder=get_contig_prefix
	log:
		"log/{sample}_eggnog.log"
	singularity:
		"docker://reslp/eggnog-mapper"
	threads: 16
	shell:
		"""
		cd results/{params.folder}
		emapper.py  -i {params.pred_folder}_preds/predict_results/{params.folder}.proteins.fa --output {params.pred_folder}_preds/eggnog_results -d euk --data_dir /data/eggnogdb --cpu {threads} --override -m diamond >& ../../{log}
		touch ../../{output}
		"""
rule annotate:
	input:
		rules.iprscan.output,
		rules.remote.output,
		rules.eggnog.output
	output:
		"results/{sample}/{sample}_annotate.done"
	params:
		folder="{sample}",
		pred_folder=get_contig_prefix
	log:
		"log/{sample}_annotate.log"
	threads: 16
	shell:
		"""
		cd results/{params.folder}
		touch ../../{output}
		funannotate annotate -i {params.pred_folder}_preds --sbt ../../data/genbank_template.txt --eggnog {params.pred_folder}_preds/eggnog_results.emapper.annotations --busco_db metazoa --cpus {threads} >& ../../{log}
		#funannotate annotate -i {params.pred_folder}_preds --sbt ../../data/genbank_template.txt --cpus {threads} >& ../../{log}
		"""
