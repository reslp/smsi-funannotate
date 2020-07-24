# A containerized genome annotation pipeline using snakemake, singularity and funannotate

snakemake -> singularity -> funannotate

This snakemake pipeline implements funannotate for use on SLURM and SGE based clusters. It was originally built to work on the VSC (Vienna Scientific Cluster with SLURM) and the University of Graz cluster (using SGE). It should be relatively simple to adopt to other cöusters (or even single servers) if you follow the guidelines below.

To make it work, a few things net to be setup beforehand:


## **Prerequisites**

- A Linux cluster or Server (Although it may also work on Apple computers)
- globally installed SLURM 18.08.7.1
- globally installed singularity 3.4.1+ 
- installed snakemake 5.15.0 (eg. in an anaconda environment)

## Rulegraph

<img src="https://github.com/reslp/smsi_funannotate/blob/master/rulegraph.png" eight="500">

## There are some issues which are specific for the respective clusters:
- Sauron (SGE): There is a strange issue with ReapeatMasker. For some reason it does not run. ReapeatMasking with TanTan works fine.
- Sauron (SGE): I have had many jobs failing due to a singularity error: `FATAL:   container creation failed: failed to resolved session directory`. This does not occur on VSC. From the extended message: `Activating singularity image /cl_tmp/reslph/projects/xylographa_fun/.snakemake/singularity/195cc8bdbe1d3f304062822f8f4f06ce.simg
FATAL:   container creation failed: failed to resolved session directory /usertmp/singularity/mnt/session: lstat /tmp/singularity: no such file or directory` I assume it has to do with the tmp directory not being present. I have seen this after the jobs have been in the queue for a week (and other jobs ran fine). Maybe the /tmp directory is automatically deleted from time to time which causes this error.

## **Pipeline setup and configuration**

### 1. Clone the repository:

The first step is to clone this repository:

```
git clone https://github.com/reslp/smsi-funannotate.git
```



### 2. Setup of funannotate database:

Funannotate requires a database containing things like BUSCO sets, PFAM information, data from dbcan etc. To setup the database, it is necessary to bind the external database directory into the container to the correct mountpoint. Run this command inside the just cloned github directory:

```
singularity shell docker://reslp/funannotate:1.7.4 -B $(pwd)/data/external:/data/database
```

You will end up inside the container. To setup and check the funannotate database run these two commands inside the container:

```
funannotate setup -i all
funannotate database
```

The database needs to be bound every time funannotate is run. This should be done automatically by the submission script (submit.sh). So that the script is able to find the database, it is important that the predefined paths are used.


### 3. Configure external dependencies:

Funannotate needs some additional programs which can not be included in the container. These include SignalP4.1, Genemark ES, interproscan and eggnog-mapper, which require you to get a license and download them manually. In the case of Eggnog mapper I have prepared a seperate container which can be used.


#### 3.1 Eggnog mapper:

Eggnog mapper comes as another container named reslp/eggnog-mapper:1.0.3. Eggnog Mapper V2 is not yet compatible with funannotate. For the container to work it is necessary to download NOG databases. To do so, these command needs to be run inside the data folder of the current project:

```
cd data
mkdir eggnogdb
singularity run docker://reslp/eggnog-mapper:1.0.3 download_eggnog_data.py NOG -y --data_dir eggnogdb
```

#### 3.2 GeneMark-ES:

GeneMark-ES can be downloaded here: [topaz.gatech.edu/GeneMark/license_download.cgi](http://topaz.gatech.edu/GeneMark/license_download.cgi)

Unzip the downloaded file into data/external/gm_et_linux_64. Place the license key file .gm_key into your project folder, at the same level where the snakefile is.

Tested with Version 4.


#### 3.3 InterProScan:

Download interproscan from [www.ebi.ac.uk/interpro/download/](https://www.ebi.ac.uk/interpro/download/). Place it in the folder data/external/interproscan-versionXXX. Make sure the path in the Snakefile points to the correct directory. InterProScan is frequently updated and your version could be different from the one specifid in the Snakefile.

IMPORTANT: Newer version of Interproscan are not compatible with funannotate <1.7.3.

Tested with Version 5.39-77.0

#### 3.4 Signal-P:

Download Signal-P from [services.healthtech.dtu.dk/service.php?SignalP-5.0](https://services.healthtech.dtu.dk/service.php?SignalP-5.0). Place it in the folder data/external/signalp-4.1. Make sure the path in the Snakefile points to the correct directory. 
You also need to change the `signalp` script to point to the correct directory. It should look like this:

```
###############################################################################
#               GENERAL SETTINGS: CUSTOMIZE TO YOUR SITE
###############################################################################

# full path to the signalp-4.1 directory on your system (mandatory)
BEGIN {
    $ENV{SIGNALP} = '/data/external/signalp-4.1';
}

# determine where to store temporary files (must be writable to all users)
my $outputDir = "/tmp";

# max number of sequences per run (any number can be handled)
my $MAX_ALLOWED_ENTRIES=100000;
```


Tested with Version 4.1

#### Optional: 3.5 Repbase Repeatmasker Library for RepeatMasker

The RepBase repeat library has become prorietory. By default the funannotate containers used in the pipeline will use the built in library shipped with RepeatMasker. It is however possible to use old versions of the RepBase library. This is handled by mounting the respective directory (Libraries directory in the RepeatMasker directory) into the container like so:

```
-B /path/RepeatMasker/Libraries/:/software/RepeatMasker/Libraries
``` 

In the context of the pipeline this needs to be handled by the submission script `submit.sh`. Look there to see how to add this correctly. I usually keep a symlink of that directory in my data folder and bindmount the symlink. This is how this looks in my data directory:
```
drwxr-xr-x  5 reslp domainusers   15 21. Jul 10:36 .
drwxr-xr-x 12 reslp domainusers   30 24. Jul 21:02 ..
drwxr-xr-x  3 reslp domainusers   96 15. Jul 11:53 assemblies
-rw-r--r--  1 reslp domainusers  645 15. Jul 12:22 cluster_config-sauron.yaml
-rw-r--r--  1 reslp domainusers  838 15. Jul 11:42 cluster_config-vsc4.yaml
-rw-r--r--  1 reslp domainusers  438 17. Apr 08:42 cluster_config.yaml
-rw-r--r--  1 reslp domainusers  783 21. Jul 10:36 config.yaml
-rw-r--r--  1 reslp domainusers 8,5K 17. Apr 08:42 data_bak.csv
drwxr-xr-x  6 reslp domainusers   40 15. Jul 11:42 database
-rw-r--r--  1 reslp domainusers 9,6K 16. Jul 09:57 data.csv
-rw-r--r--  1 reslp domainusers  395 15. Jul 11:42 data_test.csv
lrwxrwxrwx  1 reslp domainusers   34 17. Apr 09:30 eggnogdb -> /cl_tmp/reslph/databases/eggnogdb/
drwxr-xr-x  8 reslp domainusers   10 15. Jul 11:42 external
-rw-r--r--  1 reslp domainusers 1,6K 21. Jul 09:47 genbank_template.txt
lrwxrwxrwx  1 reslp domainusers   48 17. Apr 09:31 RepeatMaskerLibraries -> /cl_tmp/reslph/databases/RepeatMasker/Libraries/
```

### 4. Preparing data files and specific setting

First make sure that your data.csv file is set up correctly. Review the provided example.

Make sure to configure the cluster parameters correctly (max_walltime, qos, etc.). This can be very different from the example provided in the data directory.


Get a personalized GenBank template [https://submit.ncbi.nlm.nih.gov/genbank/template/submission/](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/) and place it in data/genbank_template.txt. 


## Run the pipeline on a cluster

Have a look at the run script:

```
$ ./submit.sh
Welcome to the pipeline submission script. A script helps to submit jobs to SLURM and SGE clusters with snakemake and singularity

Usage: ./submit.sh [-v] [-t <cluster_type>] [-c <cluster_config_file>] [-s <snakemke_args>]

Options:
  -t <cluster_type> Specify explicitly which submission system is available. Options: sge, slurm
  -c <cluster_config_file> Path to cluster config file in YAML format (mandatory). 
  -s <snakemake_args> Additional arguments passed on to the snakemake command (optional). snakemake is run with --immediate-submit -pr --notemp --latency-wait 600 --use-singularity --jobs 1001 by default.
  -i "<singularity_args>" Additional arguments passed on to singularity (optional). Singularity is run with -B /tmp:/usertmp by default.
```

A full run (incl. submission) on a SLURM cluster would look like this:

```
./submit.sh -t slurm -c data/cluster-config-vsc4.yaml
```

## Run the pipeline without SLURM on a single machine:

This should run the pipeline on a single machine without SLURM job management present.

```
snakemake --use-singularity -p --singularity-args "-B $(pwd)/data/eggnogdb:/data/eggnogdb -B $(pwd)/data/database:/data/database -B $(pwd)/data/external:/data/external -B $(pwd)/data/RepeatMaskerLibraries:/software/RepeatMasker/Libraries"
```

