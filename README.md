# LRS-Assembler

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Reproducibility / Env](https://img.shields.io/badge/Env-reproducible-blue.svg)](https://github.com/BPRC-Bioinfo/LRS-Assembler/blob/main/envs/LRS-assembler.yaml) [![Install](https://img.shields.io/badge/Install-Quick%20start-brightgreen.svg)](docs/INSTALL.md) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)]()

A Snakemake pipeline for hybrid long-read assembly and targeted annotation that produces phased haplotypes and uses them to annotate regions of interest. LRS-Assembler accepts Oxford Nanopore (ONT) and PacBio (HiFi) reads, supports multiple library types (cDNA, gDNA, protein libraries), and is designed to be reproducible and accessible to scientists.

## Why this pipeline
- Produces phased haplotypes and uses them to improve annotation of highly variable or complex loci.
- Works with different library types: cDNA, gDNA, protein.
- Can perform whole-sample assembly and targeted region extraction/annotation.

## Requirements
- Linux or macOS
- git
- conda or mamba
- snakemake
- Recommended resources: 16+ CPU cores and 64+ GB RAM for large genomes; much less for small targeted-region runs.

## Quick start
1. Clone repository:
   ```
   git clone git@github.com:BPRC-Bioinfo/LRS-Assembler.git
   cd LRS-Assembler
   ```
2. Create and activate the environment:
   ```
   conda env create -f envs/LRS-assembler.yaml
   conda activate lrs_pipe
   ```
   
3. Edit the configuration for your run:
   `configs/run-config.yaml`

4. Preview the planned workflow (dry run):
   This shows which rules will run without executing them.
   ```
   snakemake -n -s scripts/LRS-assembler.smk --use-conda --cores 4
   ```
5. Run the pipeline:
   ```
   snakemake --cores 16 -s scripts/LRS-assembler.smk --use-conda
   ```
   
## Configuration overview
Open `configs/run-config.yaml` and set sample species, paths, references, and regions of interest.

Example run-config file

```yaml
species: "Scientific species name"

annotation_ref:
  Reference_name:
    accession_number: "NCBI accession"
    genome: "/path/to/local/reference/genomic.fna"
    gff: "/path/to/local/reference/genome.gff"

scaffold_ref:
  Reference_name:
    accession_number: "NCBI accession"
    genome: "/path/to/local/reference/genomic.fa"
    chr_info: "/path/to/local/info_file"

busco: "busco database"

region:
  RegionA:
    left_flank: "NCBI gene name"
    left_flank_local: "/path/to/left.fasta"
    right_flank: "NCBI gene name"
    right_flank_local: "/path/to/right.fasta"
    cDNA_library: "path/to/cDNA_reference.fasta" 
    gDNA_library: "path/to/gDNA_reference.fasta"
    protein_library: "path/to/protein_reference.fasta"
  RegionB:
    left_flank: "NCBI gene name"
    left_flank_local: "/path/to/left.fasta"
    right_flank: "NCBI gene name"
    right_flank_local: "/path/to/right.fasta"
    cDNA_library: "path/to/cDNA_reference.fa" 
    protein_library: "path/to/protein_reference.fa"

nanopore:
  Sample_1:
    - "/path/to/local/Sample_1/nanopore/directory"
    - "/path/to/local/Sample_1/nanopore/directory2"
  Sample_2:
    - "/path/to/local/Sample_2/nanopore/directory"

pacbio:
  Sample_1:
    - "/path/to/local/Sample_1/pacbio/directory"
  Sample_2:
    - "/path/to/local/Sample_2/pacbio/directory"
```

## Key fields:

- species:
  - Scientific name of the species under study.

- annotation_ref
  - Reference used for annotation. You may provide an NCBI accession (the pipeline will attempt to download it) or local paths to the genome and annotation files.
- scaffold_ref:
  - Reference used for scaffolding contigs. You may provide an NCBI accession (the pipeline will attempt to download it) or local paths.
  - NOTE: If you use a local scaffold reference, you must supply a custom chromosome info file (tab-separated) that maps FASTA headers/accessions to chromosome names (example below).  


   ```text
   accession	chromosome
   chr1	01
   chr2	02
   chr3	03
   chrX	X
   chrY	Y
   chrMT	MT
   ```
- region:
  - Specify the region name.
  - Define regions by left and right flanking gene names (the pipeline will download sequences from NCBI based on the specified species) or by local FASTA files.
  - For each region, you may provide any combination of cDNA, gDNA, and protein libraries.

- nanopore / pacbio:
  - Per-sample lists of directories containing FASTQ (ONT) or BAM (PacBio HiFi).

## Example commands
- Dry run (preview):
  snakemake -n -s scripts/LRS-assembler.smk --use-conda --cores 4

- Small test:
  snakemake --cores 4 -s scripts/LRS-assembler.smk --use-conda

- Full run:
  snakemake --cores 32 -s scripts/LRS-assembler.smk --use-conda

### Cluster usage
Snakemake integrates with schedulers. Use --profile or --cluster as appropriate to your HPC environment. See Snakemake docs for examples.

## Outputs
For each sample run you will find a results/<sample_name>/ directory that includes:

```
.
└── results/
    └── sample_name/
        ├── annotation/
        │   └── Annotation results
        ├── benchmarks/
        │   └── Time and RAM usage per rule
        ├── hifiasm/
        │   └── All data generated by hifiasm
        ├── info/
        │   └── A HTML report including assembly information.
        ├── midway/
        │   └── Intermediate working files
        ├── logs/
        │   └── Log run for some rules
        ├── raws/
        │   └── Size-filtered raw HiFi and ONT reads
        └── scaffolds/
            └── All final and intermediate scaffolds generated by LRS-Assembler. 
```

# Standalone annotation tool

A major feature of this project is a standalone annotation workflow in the annotation/ directory. This makes the annotation step usable independently (for example, if you already have assemblies).

To use:
- See the README in the annotation/ directory.
- Prepare configs/config.yaml (a template is included).
- Place assembly FASTA(s) in `annotation/inputs/` named like: {sample}_{hap}.fa
- Run:

```
  conda activate lrs_pipe
  cd annotation
  snakemake --cores 4 --use-conda
```

Troubleshooting & support
- Check results/<sample>/logs/ for failing rules.
- Re-run failed/incomplete rules:
  ```
  snakemake --cores 8 -s scripts/LRS-assembler.smk --use-conda --rerun-incomplete
  ```
  
- For quick diagnosis, use the dry run (-n) to see planned steps.
- Contact: le@bprc.nl or bruijnesteijn@bprc.nl
- To report bugs: open an issue and include the minimal config and logs (remove private data).

Cite
Please cite the project and include the version or commit hash used. Manuscript describing LRS-Assembler is in preparation.
