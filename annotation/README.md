# LRS-Assembler - Annotation 

The annotation tool can be run independently of the assembly pipeline. This process requires an assembled genome or haplotype sequence, which will first be scanned for the presence of flanking genes. Once identified, the annotation process begins using the user-provided reference library.

## Installation Guide

To download the annotation tool, the repository can be cloned using the following command:

```
git clone git@github.com:BPRC-Bioinfo/LRS-Assembler.git
```
And then navigate to the stand-alone annotation tool:
```
cd LRS-Assembler/annotation

```

## Usage Guide

Like the assembly pipeline, the annotation tool also requires a configuration file. A template for this file is available: ```configs/anno_run-config.yaml```

```
species: "scientific name species"
region:
  Region1:
    left_flank: "flanking gene"
    right_flank: "flanking gene"
    library: "/path/to/references.fasta" 
    minimap2: "-cx splice:hq -G16k"
    blast: "-word_size 7"
  Region2:
    left_flank: "flanking gene"
    right_flank: "flanking gene"
    library: "/path/to/references2.fasta"
    minimap2: "-cx splice:hq -G16k"
    blast: "-word_size 7"
```

### Regions of Interest

Specify your region of interest under ```region:```.

Provide the left and right flanking genes.
If only one flanking gene is specified, the program will process the region from that gene to the end of the sequence.

### Library

The library file requires a specific format to annotate the RIO, which includes the gene name and the length of the sequence.
The expected format:

```
>Gene|length
```

```
# create environment for renaming your library
conda env create -f ../envs/rename.yaml

# Activate environment
conda activate lib_prepare

# Rename the fasta file
python ../scripts/lib_format.py input_ibrary.fa output_library.fa

# Validate the renamed library
grep ">" output_library.fa
```

Specify the path to the reference library in the configuration file.

## Prepare Input Files

Create a directory named ```inputs``` and place your assembly or sequence files inside it.

## Run the Annotation Tool

```
conda activate snakemake

snakemake --use-conda 
```
