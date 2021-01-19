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
		if [[ "{params.premasked}" == "yes" ]]; then # yes means assembly is premasked
			cp {input.assembly} {output}
		elif [[ "{params.premasked}" == "tantan" ]]; then # means that tantan should be used to mask the genome
			cd results/{params.folder}
			funannotate mask -i ../../{input.assembly} -o ../../{output} -m tantan --cpus {threads}
		else # use repeatmasker or rather the method specified in config file to mask
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
		# this is added to prevent tRNAscan from failing
		export TMPDIR="$(pwd)/tmp"

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
	singularity: "docker://reslp/interproscan-wrapper:5.48-83.0"
	threads: 16
	shell:
		"""
		cd results/{params.folder}
		mkdir -p {params.pred_folder}_preds/annotate_misc
		#funannotate iprscan --iprscan_path /data/external/interproscan-5.33-72.0/interproscan.sh -i ../../results/{params.folder}/{params.pred_folder}_preds -m local -c 2 >& ../../{log}
		/data/external/interproscan-5.48-83.0/interproscan.sh -cpu {threads} -T $(pwd)/interpro_tmp -i ../../results/{params.folder}/{params.pred_folder}_preds/predict_results/{params.folder}.proteins.fa -o ../../results/{params.folder}/{params.pred_folder}_preds/annotate_misc/iprscan.xml -f XML -goterms -pa >& ../../{log}
		rm -rf $(pwd)/interpro_tmp
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
		pred_folder=get_contig_prefix,
		buscodb=config["annotate"]["buscodb"]
	log:
		"log/{sample}_annotate.log"
	threads: config["annotate"]["threads"]
	shell:
		"""
		cd results/{params.folder}
		touch ../../{output}
		funannotate annotate -i {params.pred_folder}_preds --sbt ../../data/genbank_template.txt --eggnog {params.pred_folder}_preds/eggnog_results.emapper.annotations --busco_db {params.buscodb} --cpus {threads} >& ../../{log}
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
	
