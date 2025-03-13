# Example run with LRS-Assembler


This tutorial describes how to modify your configuration file to run LRS-Assembler on multiple samples.
If you are processing several samples with the same settings, it is best to create a dedicated configuration file.
In this example, we show how to assemble multiple rhesus macaque genomes and annotate the LILR and KIR regions.

## 1. Create a new config file


Start by copying the existing configuration file to create a new one:


```
cp configs/run-config.yaml configs/macaque.yaml
```


Next, open the new configuration file in your preferred text editor:


```
nano configs/macaque.yaml
```


If you do not have the reference files available locally, the pipeline will automatically download and format them using the provided NCBI accession number.

Below is an example configuration file:


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
    blast: "-word_size 7" # BLAST setting
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

*Remember to update the file paths to match your local directory structure.*

### Editing the Snakemake File


Update the Snakemake file so that it loads the correct configuration file.
Open the file using your favorite editor:


```
nano scripts/Snakefile
```


Ensure the following line is pointing to your new configuration file:


```
configfile: "configs/macaque.yaml"
```


### Run the pipeline:


To run the pipeline, use the command below. Replace `<number_of_cores>` with the number of CPU cores you want to use:


```
snakemake --cores <number_of_cores> -s scripts/Snakefile --use-conda
```

## 2. Custom info file


The custom info file improves liftoff’s performance by linking scaffold names in your reference to the correct chromosomes. 
The `@` symbol acts as a placeholder, allowing the pipeline to dynamically insert the sample name.
An example file is provided below.


```
NC_041754.1,@_chr01
NC_041755.1,@_chr02
NC_041756.1,@_chr03
NC_041757.1,@_chr04
NC_041758.1,@_chr05
NC_041759.1,@_chr06
NC_041760.1,@_chr07
NC_041761.1,@_chr08
NC_041762.1,@_chr09
NC_041763.1,@_chr10
NC_041764.1,@_chr11
NC_041765.1,@_chr12
NC_041766.1,@_chr13
NC_041767.1,@_chr14
NC_041768.1,@_chr15
NC_041769.1,@_chr16
NC_041770.1,@_chr17
NC_041771.1,@_chr18
NC_041772.1,@_chr19
NC_041773.1,@_chr20
NC_041774.1,@_chrX
NC_027914.1,@_chrY
```

## 3. BUSCO evaluation


Ensure you choose the appropriate BUSCO lineage database for the species you are evaluating.
For the rhesus macaque, the configuration should include:

```
busco: "primates_odb10" 
```

For more information on BUSCO databases and to determine which lineage database is best for your species, please visit the [BUSCO website](https://busco.ezlab.org).
