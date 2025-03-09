## Setup config file to run multiple samples


This tutorial explains how to modify your configuration file to run LRS-Assembler with multiple samples.
If you are running several samples with the same settings, it is advisable to create a dedicated configuration file for them.
Below is an example for assembling multiple rhesus macaque genomes and annotating the LILR and KIR regions.

### Create a new config file

First, copy the existing configuration file to create a new one:

```
cp configs/run-config.yaml configs/macaque.yaml
```

Next, edit the new configuration file using your preferred text editor:


```
species: "macaca mulatta"
reference:
  Mmul10: # Name of the reference
    accession_number: "GCA_003339765.3" # NCBI accession number of the reference
    genome: "/path/to/local/reference/genomic.fna" # local location of the genome
    gff: "/path/to/local.gff" # local location of the gff
    chr_info: "/path/to/local/chr_info_file" # local location of the chromosome
region:
  LILR: # Name of the region
    left_flank: "RPS9" # The flanking gene name
    right_flank: "FCAR" # The flanking gene name
    library: "/path/to/LILR.fa" # The reference fasta
    minimap2: "-x splice:hq" # minimap2 setting for cDNA
    blast: "-word_size 7" # blast setting
  KIR:
    left_flank: "FCAR"
    right_flank: "LILRA6"
    library: "/path/to/referenceB.fasta" 
    minimap2: "-x asm5" # minimap2 setting for gDNA
    blast: "-word_size 7"
busco: "primates_odb10" # BUSCO lineage database setting (-l)
nanopore:
  AnimalA1: # Name of the sample
    - "/path/to/local/Sample_1/nanopore/fastq1" # Directory for the Nanopore raw files (fastq.gz)
    - "/path/to/local/Sample_1/nanopore/fastq2" # Directory 2 for the Nanopore raw files (fastq.gz)
  AnimalB2:
    - "/path/to/local/Sample_2/nanopore/fastq1" 
pacbio:
  AnimalA:
    - "/path/to/local/Sample_1/pacbio/bam1" # Directory for the Pacbio raw files (.bam & .bai)
    - "/path/to/local/Sample_1/pacbio/bam2" # Directory 2 for the Pacbio raw files (.bam & .bai)
  AnimalB2:
    - "/path/to/local/Sample_2/pacbio/bam1"
```

### Editing the Snakemake File


Edit the `scripts/Snakemake` file to load in the correct config file.

```
nano scripts/Snakefile
```

Ensure the following line is pointing to your new configuration file:

```
configfile: "configs/macaque.yaml"
```

### Run the pipeline:

To run the pipeline, use the following command:


```
snakemake --cores <number_of_cores> -s scripts/Snakefile --use-conda
```


