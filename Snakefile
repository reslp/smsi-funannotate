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

rule clean:
	input:
		assembly = get_assembly_path
	
	output:
		"results/{sample}/{sample}_cleaned.fas"
	log:
		"log/{sample}_clean.log"
	params:
		folder = "{sample}",
		minlen = config["clean"]["minlen"]
	shell:
		"""
		if [[ ! -d results/{params.folder} ]]
		then
			mkdir results/{params.folder}
		fi
		cd results/{params.folder}
		funannotate clean -i ../../{input.assembly} -o ../../{output} --minlen {params.minlen}  2> ../../{log}
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
		premasked = get_premasked_state,
		method = config["mask"]["method"]
	threads: config["mask"]["threads"]
	shell:
		"""
		if [[ "{params.premasked}" == "yes" ]]; then
			cp {input.assembly} {output}
		else
			cd results/{params.folder}
			funannotate mask -i ../../{input.assembly} -o ../../{output} -m {params.method} --cpus {threads}
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
		organism = config["predict"]["organism"],
		busco_seed_species = config["predict"]["busco_seed_species"],
		ploidy = config["predict"]["ploidy"],
		busco_db = config["predict"]["busco_db"]
	log:
		"log/{sample}_predict.log"
	threads: config["predict"]["threads"] 
	shell:
		"""
		if [[ ! -d results/{params.folder} ]]
		then
			mkdir results/{params.folder}
		fi
		cd results/{params.folder}
		funannotate predict -i ../../{input.assembly} -o {params.pred_folder}_preds -s {params.sample_name} --name {params.pred_folder}_pred --optimize_augustus --cpus {threads} --busco_db {params.busco_db} --organism {params.organism} --busco_seed_species {params.busco_seed_species} --ploidy {params.ploidy} >& ../../{log}
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
		#funannotate iprscan --iprscan_path /data/external/interproscan-5.33-72.0/interproscan.sh -i ../../results/{params.folder}/{params.pred_folder}_preds -m local -c 2 >& ../../{log}
		/data/external/interproscan-5.39-77.0/interproscan.sh -i ../../results/{params.folder}/{params.pred_folder}_preds/predict_results/{params.folder}.proteins.fa -o ../../results/{params.folder}/{params.pred_folder}_preds/annotate_misc/iprscan.xml -f XML -goterms -pa >& ../../{log}
		touch ../../{output}
		"""
rule remote:
	input:
		rules.predict.output
	output:
		"results/{sample}/{sample}_remote.done"
	params:
		folder="{sample}",
		pred_folder=get_contig_prefix,
		methods = config["remote"]["methods"],
		email = config["remote"]["email"]
	log:
		"log/{sample}_remote.log"
	shell:
		"""
		cd results/{params.folder}
		funannotate remote -i {params.pred_folder}_preds -m {params.methods} -e {params.email} >& ../../{log}
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
		"docker://reslp/eggnog-mapper:1.0.3"
	threads: config["eggnog"]["threads"]
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
	threads: config["annotate"]["threads"]
	shell:
		"""
		cd results/{params.folder}
		touch ../../{output}
		funannotate annotate -i {params.pred_folder}_preds --sbt ../../data/genbank_template.txt --eggnog {params.pred_folder}_preds/eggnog_results.emapper.annotations --busco_db metazoa --cpus {threads} >& ../../{log}
		#funannotate annotate -i {params.pred_folder}_preds --sbt ../../data/genbank_template.txt --cpus {threads} >& ../../{log}
		"""

species_names, preddirs = glob_wildcards("results/{sample}/{preddir}_preds")
# funannotate compare does not allow / charcters in the output folder
# therefore the folder has to be moved manually to the results folder.
if config["compare"]["phylogeny"] == "yes" or config["compare"]["histograms"] == "yes":
	rule compare:
                input:
                        checkpoint=expand("results/{sam}/{sam}_annotate.done", sam=sample_data.index.tolist()),
                        folders = expand("results/{species_name}/{preddir}_preds", zip, species_name=species_names, preddir=preddirs)
                output:
                        checkpoint = "results/funannotate_compare.done",
                        dir = directory("results/funannotate_compare/")
                params:
                        samples = expand("results/{sam}", sam=sample_data["contig_prefix"].tolist()),
                        num_orthos = config["compare"]["num_orthos"],
                        ml_method = config["compare"]["ml_method"]
                singularity:
                        "docker://reslp/funannotate:1.7.2"
                log:
                        "log/funannotate_compare.log"
                threads: config["compare"]["threads"]
                shell:
                        """
                        funannotate compare --cpus {threads} --num_orthos {params.num_orthos} --ml_method {params.ml_method} -i {input.folders} >& {log}
                        cp -r funannotate_compare results/
                        cp funannotate_compare.tar.gz results/
                        rm -rf funannotate_compare
			rm funannotate_compare.tar.gz
			touch {output.checkpoint}
                        """	
else:
	rule compare:
		input:
			checkpoint=expand("results/{sam}/{sam}_annotate.done", sam=sample_data.index.tolist()),
			folders = expand("results/{species_name}/{preddir}_preds", zip, species_name=species_names, preddir=preddirs) 
		output:
			checkpoint = "results/funannotate_compare.done",
			dir = directory("results/funannotate_compare/")
		params:
			samples = expand("results/{sam}", sam=sample_data["contig_prefix"].tolist()),
			num_orthos = config["compare"]["num_orthos"],
			ml_method = config["compare"]["ml_method"]
		singularity:
			"docker://reslp/funannotate:experimental"
		log:
			"log/funannotate_compare.log"
		threads: config["compare"]["threads"]
		shell:
			"""
			funannotate compare --cpus {threads} --num_orthos {params.num_orthos} --ml_method {params.ml_method} -i {input.folders} >& {log}
			cp -r funannotate_compare results/
                        cp funannotate_compare.tar.gz results/
                        rm -rf funannotate_compare
                        rm funannotate_compare.tar.gz
			touch {output.checkpoint}
			"""

def get_ids_for_idsfile(wildcards):
	prefix = sample_data["contig_prefix"].tolist()
        species = sample_data["sample"].tolist()
	ids = ""
	for pre, sp in zip(prefix,species):
		ids = ids + pre +"\t"+ sp + "\n" 
	ids = ids.strip("\n")
	return ids

rule prepare_downstream_output:
		input:
			gff_files = expand("results/{sam}/{pref}_preds/annotate_results/{sam}.gff3", zip, sam=species_names, pref=preddirs),
			protein_files = expand("results/{sam}/{pref}_preds/annotate_results/{sam}.proteins.fa", zip, sam=species_names, pref=preddirs),
			transcript_files = expand("results/{sam}/{pref}_preds/annotate_results/{sam}.cds-transcripts.fa", zip, sam=species_names, pref=preddirs),
			cazy_results = "results/funannotate_compare/cazy/CAZyme.all.results.csv",
			cazy_summary_results = "results/funannotate_compare/cazy/CAZyme.summary.results.csv",
			interproscan_results = "results/funannotate_compare/interpro/interproscan.results.csv",
			pfam_results = "results/funannotate_compare/pfam/pfam.results.csv"
		output:
			combined_gff = "results/downstream/combined.gff",
			combined_proteins = "results/downstream/combined_proteins.fa",
			combined_transcripts = "results/downstream/combined_transcripts.fa",
			cazy_results = "results/downstream/CAZyme.all.results.csv",
			cazy_summary_results = "results/downstream/CAZyme.summary.results.csv",
			interproscan_results = "results/downstream/interproscan.results.csv",
			pfam_results = "results/downstream/pfam.results.csv",
			ids = "results/downstream/ids.txt",
			protein_files = directory("results/downstream/protein_files/"),
			gff_files = directory("results/downstream/gff_files/")
		params:
			ids = get_ids_for_idsfile
		shell:
			"""
			echo "##gff-version 3" > {output.combined_gff}
			sed '/##gff-version 3/d'  {input.gff_files} >> {output.combined_gff}
			cat {input.protein_files} > {output.combined_proteins}
			cat {input.transcript_files} > {output.combined_transcripts}
			cp {input.cazy_results} {output.cazy_results}
			cp {input.cazy_summary_results} {output.cazy_summary_results}
			cp {input.interproscan_results} {output.interproscan_results}
			cp {input.pfam_results} {output.pfam_results}
			mkdir -p {output.gff_files}
			cp {input.gff_files} {output.gff_files}
			mkdir -p {output.protein_files}
			cp {input.protein_files} {output.protein_files}
			echo "{params.ids}" > {output.ids}
			"""
	
