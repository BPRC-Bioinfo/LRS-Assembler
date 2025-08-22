# LRS-Assembler

[![CI](https://github.com/BPRC-Bioinfo/LRS-Assembler/actions/workflows/ci.yml/badge.svg)](https://github.com/BPRC-Bioinfo/LRS-Assembler/actions) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Reproducibility / Env](https://img.shields.io/badge/Env-reproducible-blue.svg)](https://github.com/BPRC-Bioinfo/LRS-Assembler/blob/main/envs/LRS-assembler.yaml) [![Install](https://img.shields.io/badge/Install-Quick%20start-brightgreen.svg)](docs/INSTALL.md) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)]()

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
   git clone git@github.com:BPRC-Bioinfo/LRS-Assembler.git
   cd LRS-Assembler

2. Create and activate the environment:
   conda env create -f envs/LRS-assembler.yaml
   conda activate lrs_pipe
   (Or use mamba for faster installs:
   mamba env create -f envs/LRS-assembler.yaml)

3. Edit the configuration for your run:
   configs/run-config.yaml

4. Preview the planned workflow (dry run):
   snakemake -n -s scripts/LRS-assembler.smk --use-conda --cores 4
   This shows which rules will run without executing them.

5. Run the pipeline:
   snakemake --cores 16 -s scripts/LRS-assembler.smk --use-conda

## Configuration overview
Open configs/run-config.yaml and set sample species, paths, references, and regions of interest.

Example minimal config snippet
```yaml
species: "Homo sapiens"

annotation_ref:
  GRCh38:
    accession_number: "GCF_000001405.39"
    genome: "/path/to/GRCh38.fa"
    gff: "/path/to/GRCh38.gff"

region:
  MHC:
    left_flank_local: "/path/to/HLA-A.fa"
    right_flank_local: "/path/to/HLA-B.fa"
    cDNA_library: "/path/to/cDNA_refs.fa"

nanopore:
  sample_1:
    - "/data/sample_1/nanopore/"

pacbio:
  sample_1:
    - "/data/sample_1/pacbio/"
```

## Key fields:

- species:
  - Scientific name of the species under study.

- annotation_ref / scaffold_ref:
  - You may provide an NCBI accession (the pipeline will try to download) or local paths to genome and annotation files.
  - NOTE: If you choose to use a local scaffold reference, you must supply a custom chromosome info file (tab-separated) that maps FASTA headers/accessions to chromosome names (example below).

accession<TAB>chromosome
chr1<TAB>01
chr2<TAB>02
chr3<TAB>03
chrX<TAB>X
chrY<TAB>Y
chrMT<TAB>MT

(Example lines as a file)
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
  - Define regions by left/right flanking gene names or by local FASTA files for the flanks.
  - For each region, you can provide cDNA, gDNA, and protein libraries.

- nanopore / pacbio:
  - Per-sample lists of directories containing FASTQ/GZ (ONT) or BAM/FASTQ (PacBio HiFi).

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
- annotation/ — final annotation results
- benchmarks/ — time & RAM usage
- hifiasm/ — hifiasm outputs
- info/ — HTML report with assembly stats
- midway/ — intermediate files
- logs/ — program logs
- raws/ — size-filtered reads
- scaffolds/ — final and intermediate scaffolds

# Standalone annotation tool

A major feature of this project is a standalone annotation workflow in the annotation/ directory. This makes the annotation step usable independently (for example, if you already have assemblies).

To use:
- Prepare configs/config.yaml (template included).
- Place assembly FASTA(s) in annotation/inputs/ named like: {sample}_{hap}.fa
- Run:

```
  conda activate lrs_pipe
  cd annotation
  snakemake --cores 4 --use-conda
```

Troubleshooting & support
- Check results/<sample>/logs/ for failing rules.
- Re-run failed/incomplete rules:
  snakemake --cores 8 -s scripts/LRS-assembler.smk --use-conda --rerun-incomplete
- For quick diagnosis, use the dry run (-n) to see planned steps.
- Contact: le@bprc.nl or bruijnesteijn@bprc.nl
- To report bugs: open an issue and include the minimal config and logs (remove private data).

Cite
Please cite the project and include the version or commit hash used. Manuscript describing LRS-Assembler is in preparation.

License
Add or confirm a LICENSE file (recommended: MIT). If you tell me which license you prefer, I can add the LICENSE file.
