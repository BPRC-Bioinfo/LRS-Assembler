# Instruction

Download the git repo.

```
git clone git@github.com:BPRC-Bioinfo/LRS-Assembler.git

## Go to the standalone annotation directory
cd LRS-Assembler/annotation

```

## Modify the config:

A template of the config file is in configs/anno_run-config.yaml

```
species: "scientific species name"
region:
  Region1:
    left_flank: "flanking gene"
    right_flank: "flanking gene"
    library: "/path/to/references.fasta" 
    minimap2: "-cx splice:hq -G16k"
    blast: "-word_size 7"
  Region2:
    left_flank: "flanking gene"
    right_flank: ""
    library: "/path/to/references2.fasta"
    minimap2: "-cx splice:hq -G16k"
    blast: "-word_size 7"
```

Type in the scientific species name.

### Regions of interest

Name for your region of interest underneath `regions:`


Fill in the left and right flanking genes.
If only one flanking gene is filled, the program will process the region from the flanking gene to the end of the sequence.

### Library

The program needs the length of the fasta sequence so that it can do additional calculation.
This is what it should look like

```
>Gene|length
```

```
# create environment for renaming
conda env create -f ../envs/rename.yaml

# Activate environment
conda activate lib_prepare

# Rename the fasta file
python ../scripts/lib_format.py input_ibrary.fa output_library.fa

# Check if it works
grep ">" output_library.fa
```

Add path of the reference library to the config file. 

## Preparing analysis sequences

Create a directory called `inputs`
Plase your `fasta_sequence.fa` files in the directory `inputs`

## Running the annalysis

```
conda activate snakemake

snakemake --use-conda 
```
