singularity: "docker://reslp/funannotate:1.7.2"

configfile: "data/config.yaml"

import pandas as pd
import os

include: "rules/utilities.smk"

rule all:
	input:
		#expand("results/{name}/{name}_remote.done", name=sample_data.index.tolist()),
		#expand("results/{name}/{name}_iprscan.done", name=sample_data.index.tolist()),
		#expand("results/{name}/{name}_eggnog.done", name=sample_data.index.tolist()),
		expand("results/{name}/{name}_tarpredict.done", name=sample_data.index.tolist()),
		expand("results/{name}/{name}_predict.done", name=sample_data.index.tolist()),
		expand("results/{name}/{name}_annotate.done", name=sample_data.index.tolist()),
		"results/funannotate_compare.done"

include: "rules/prepare_assembly.smk"
include: "rules/funannotate.smk"
