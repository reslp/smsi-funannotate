# funannotate -> singularity -> snakemake

### Used snakemake command:
```
snakemake --use-singularity --jobs 100 --cluster-config data/cluster_config.yaml --cluster '/home/lv71312/reslp/projects/monos_fun/bin/immediate_submit.py {dependencies}' --immediate-submit -pr --notemp
```

