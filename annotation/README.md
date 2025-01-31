

# Instruction

Download the branch annotation

```
git clone -b annotation git@github.com:BPRC-Bioinfo/LRS-Assembler.git

cd LRS-Assembler/annotation

```

Modify the config:

```
species: "macaca mulatta"
region:
  KIR:
    left_flank: "FCAR"
    right_flank: "LILRA6"
    library: "references/mamu_kir_gen_2501.fasta" 
    minimap2: "-cx splice:hq -G16k"
    blast: "-word_size 7"
  LILR:
    left_flank: "RPS9"
    right_flank: "FCAR"
    library: "references/LILR_cDNA_library.fasta"
    minimap2: "-cx splice:hq -G16k"
    blast: "-word_size 7"
```

Type in the species.
Fill in the flanking genes
Add path of the reference file


Your input files, should be in the directory `inputs`

to run the annalysis

```
conda activate snakemake

snakemake -s snakefile.py --use-conda 


```
