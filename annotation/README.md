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
    left_flank: "Flanking gene"
    left_flank_local: "/path/to/left.fasta"
    right_flank: "Flanking gene"
    right_flank_local: "/path/to/right.fasta"
    cDNA_library: "/path/to/cDNA_reference.fasta" 
    gDNA_library: "/path/to/gDNA_references.fasta"

  Region2:
    left_flank: "Flanking gene"
    left_flank_local: "/path/to/left.fasta"
    right_flank: "Flanking gene"
    right_flank_local: "/path/to/right.fasta"
    cDNA_library: "/path/to/cDNA_reference.fasta"
```

### Regions of Interest

Specify the name of your region of interest under ```region:```.

Provide the left and right flanking genes.
If only one flanking gene is specified, the program will process the region from that gene to the end of the sequence.

### Library

Specify the path of the reference library in the configuration file.
You can specify cDNA and/or gDNA libraries.
Only list the library that you have.
Duplicate records in the library will be removed.

## Prepare Input Files

Create a directory named ```inputs``` and place your assembly or sequence files inside it.
The program can detect and process `.fa` or `.fasta` files.

## Run the Annotation Tool

```
conda activate snakemake

snakemake --use-conda 
```

The final results can be found inside the `LRS-annotation` directory.
