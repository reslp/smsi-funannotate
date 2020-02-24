# A containerized snakemake funannotate pipeline 

funannotate -> snakemake -> singularity


This snakemake pipeline implements funannotate for use on SLURM based clusters. Currently it is set up to work on the VSC (Vienna Scientific Cluster). It should be relatively simple to adopt it on other SLURM clusters or even on cluster which use a different jobs scheduling system such as SGE.

To make it work, a few things net to be setup beforehand:


## **Prerequisites**:

- A Linux cluster
- globally installed SLURM 18.08.7.1
- globally installed singularity 3.4.1+ 
- installed snakemake 5.10.0 (can be in an anaconda environment)




## **Setup of funannotate database:**

For the database setup, it is necessary to bind the external database directory into the container to the correct mountpint. I singularity like this:

	singularity shell -B /external/path:/data/database

Inside the container the funannotate db can now be set up:

	funannotate setup -i all
	funannotate database

The database needs to be bound every time funannotate is run. This should be done automatically by the submission script (bin/immediate_submitt.py).


## **External dependencies for funannotate:**

These include the programs needed for funannotate to run but which are not included in the container. It currently includes SignalP4.1, Genemark ES, interproscan and eggnog-mapper.


### Eggnog mapper:

Eggnog mapper  comes as another container named reslp/eggnog-mapper.  For it to work it is necessary to download NOG databases. To do so, these command needs to be run inside the data folder of the current project:

	mkdir eggnogdb
	singularity run docker://reslp/eggnog-mapper download_eggnog_data.py NOG -y --data_dir eggnogdb


### GeneMark-ES:

GeneMark-ES can be downloaded here: [topaz.gatech.edu/GeneMark/license_download.cgi](http://topaz.gatech.edu/GeneMark/license_download.cgi)

Unzip the downloaded file into data/external/gm_et_linux_64. Place the license key file .gm_key into your project folder, at the same level where the snakefile is.

Tested with Version 4.


### InterProScan:

Download interproscan from [www.ebi.ac.uk/interpro/download/](https://www.ebi.ac.uk/interpro/download/). Place it in the folder data/external/interproscan-versionXXX. Make sure the path in the Snakefile points to the correct directory. InterProScan is frequently updated and your version could be different from the one specifid in the Snakefile.

Tested with Version 5.39-77.0


### Signal-P:

Download Signal-P from [services.healthtech.dtu.dk/service.php?SignalP-5.0](https://services.healthtech.dtu.dk/service.php?SignalP-5.0). Place it in the folder data/external/signalp-4.1. Make sure the path in the Snakefile points to the correct directory. 

Tested with Version 4.1


## **Preparing data files and specific setting**

First make sure that your data.csv file is set up correctly. Review the provided example.

Make sure to configure the cluster parameters currectly (max_walltime, qos, etc.). This can be very different from the example provided in the data directory.


Get a personalized GenBank template [https://submit.ncbi.nlm.nih.gov/genbank/template/submission/](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/) and place it in data/genbank_template.txt. 


## **Run the pipeline**

A dry run can be started with:

	./submit.sh -n

A full run (incl. submission) like this:

	./submit.sh
