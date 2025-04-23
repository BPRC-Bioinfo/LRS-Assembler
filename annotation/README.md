# LRS-assembler - Standalone annotation tool 

The annotation tool can be run independently of the assembly pipeline. 
This process requires assembled genome or haplotype sequences, which will first be scanned for flanking genes. 
Once identified, the annotation proceeds using the user-provided reference library.

## Installation

To download the annotation tool, clone the repository using the following command:

```
git clone git@github.com:BPRC-Bioinfo/LRS-Assembler.git
```

Navigate to the stand-alone annotation tool directory:

```
cd LRS-Assembler/annotation
```

## Usage

The annotation tool, like the assembly pipeline, requires a configuration file. A template for this file is available at: ```configs/anno_run-config.yaml```

```
species: "scientific name of species"
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


### Regions of interest

Rename ```Region1``` to the specific name of your region of interest.
You can define multiple regions, each with distinct flanking genes and reference libraries.


### Flanking genes

Provide the names of the left and right flanking genes.
The program will attempt to download these genes from NCBI for the specified scientific species. Alternatively, you can provide a local FASTA file for each flanking gene if preferred (in `.fasta` format).

### Library

Specify the file paths for your reference cDNA and/or gDNA libraries in the configuration file.
Only include the paths for the libraries you intend to use.
Duplicate records within the provided libraries will be automatically removed.

## Prepare Input Files

Create a directory named ```inputs``` and place your assembly or sequence files (in `.fa` or `.fasta` format) within it.
Ensure that your sample files are named according to the following format:

```
{sample}_{hap}.fa or {sample}_{hap}.fasta
```

Where ```{sample}``` represents the name of your sample, and ```{hap}``` indicates the haplotype (e.g., hap1, hap2).
This naming convention allows the program to group haplotypes belonging to the same sample in the final report.


## Run the Annotation Tool

```
conda activate snakemake

snakemake -c {core} --use-conda 
```

The final results can be found inside the `LRS-annotation` directory, which will contain a html report file.
