# Calling HPV Integration Using Illumina Short-Read Sequencing
Workflow for calling HPV integration sites and events in Illumina short-read sequencing data. 

# Installation
This will clone the repository. You can run workflow within this directory.
```
git clone https://github.com/vanessa-porter/callRNAHPVFusions
```

### Dependencies
> To run this workflow, you must have snakemake (v6.12.3) and conda. You can install snakemake using [this guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). The remaining dependencies will be downloaded automatically within the snakemake workflow.

# Input Files

### RNA-seq <br />
- RNA alignment (bam file)
- HPV integration sites (txt file)

# Set Up Configuration Files

### **Edit the config files**

#### **Example parameters.yaml:** <br />
Config files to specify parameters and paths needed for the workflow. The main parameter to include is the genome path and the STAR transcriptome reference directory.

```
genome_path: /path/to/genome/fasta
transcriptome_path: "/path/to/star/transcriptome/directory"
```

#### **Example samples.yaml:** <br />
Main config file to specify input files.

```
samples:
    sampleName_1: 
        ont-dna-sites: /path/to/hpv_integration_sites.bed
        rna: /path/to/rnaseq.bam
    sampleName_2: 
        ont-dna-sites: /path/to/hpv_integration_sites.bed
        rna: /path/to/rnaseq.bam
```

#### Converting sample paths to yaml file
A text file can be converted to the samples.yaml file using the scripts/samplestoyaml.py script. The tsv file should be tab delimited and have the sample name in the first column, the file type in the second column (ont-dna-sites, rna) and the path in the third column (no header). 

```
scripts/sampletsvtoyaml.py -t samples.txt -o config/samples.yaml
```

# **Run Workflow**
This is the command to run the workflow with snakemake using conda. The `-c` parameter can be used to specify maximum number of threads. 

```
snakemake -c 30 --use-conda
```

# Contributors
The contributors of this project are Vanessa Porter

<a href="https://github.com/vanessa-porter/illuminaCallHPVInt/graphs/contributors">
</a>