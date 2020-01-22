configfile: "data/config.yaml"

import pandas as pd

samples = pd.read_table("data/data.csv").set_index("sample", drop=False)
print(samples.assembly_path.tolist())


rule all:
	input:
		expand("results/output_{name}.txt", name=samples.index.tolist())

rule clean:
	output:
		"results/output_{sample}.txt"
	shell:
		"""
		touch {output}
		"""