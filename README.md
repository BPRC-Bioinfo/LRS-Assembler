# Long Read Sequencing - Assembler

This Snakemake workflow assembles long reads from Oxford Nanopore (ONT) and PacBio (HiFi) datasets to create hybrid assemblies.
The workflow is designed for scalability and flexibility, making it ideal for genome assembly projects where both ONT and PacBio data are available.

## Features

  * Support for ONT and PacBio: Assemble genomes using both sequencing platforms for a hybrid assembly.
  * Automatic Sample Detection: The pipeline automatically detects input data, ensuring it only processes samples with both ONT and PacBio reads.
  * Configurable Assembly Options: Users can customize the workflow by providing references for alignment and comparison.
  * Reference Comparison: Align assembled sequences to reference genomes and compute assembly statistics.

## Workflow Overview

(Figure)

The workflow automatically detects samples and processes only those that contain both ONT and PacBio data, ensuring hybrid assembly:

  * Input Files Detection: The workflow searches for ONT .fastq.gz and PacBio .bam files in the designated directories.
  * Hybrid Assembly: The pipeline performs assembly using tools optimized for combining ONT and PacBio reads for higher assembly accuracy.
  * Contig lengthening using long reads
  * Scaffolding based on the reference genome provided
  * Evaluate using BUSCO and Quast
  * Post-assembly Processing: The assembled sequences are aligned to multiple references and statistics are computed for evaluation.

Requirements

    Snakemake (>=6.0)
    Pandas 2.2.3
    Python 3.12.7

Installation

    Clone this repository:

    git clone git@github.com:BPRC-Bioinfo/LRS-Assembler.git
    cd LRS-Assembler

    Install environment:
    ```
    conda env create -f envs/LRS-assembler.yaml
    conda activate lrs_pipe
    ```
    Customize the configs/run-config.yaml file with your reference genome information and other settings.

Input Data Structure

The input data should be organized as follows:

```
raws/
├── sample1/
│   ├── nanopore/
│   │   └── reads.fastq.gz
│   └── pacbio/
│       └── reads.bam
├── sample2/
│   ├── nanopore/
│   │   └── reads.fastq.gz
│   └── pacbio/
│       └── reads.bam
```
Configuration

Here is an example run-config.yaml file:

```
species: "macaca mulatta"
region:
  KIR:
    left_flank: "FCAR"
    right_flank: "LILRA6"
    library: "references/mamu_kir_gen_2501.fasta" 
    minimap2: "-cx splice:hq -G16k"
    blast: "-word_size 7"
```

Edit the configs/run-config.yaml file to specify:

    Library: Paths to reference files used for alignment.
    Minimap2: Customize command for minimap2
    blast: Customize command for blastn

Execution

Run the workflow with Snakemake:

```snakemake --cores <number_of_cores> -s scripts/Snakefile```


Output

The output files are stored in the results/ directory. For each sample and assembly method, you will obtain:

    Assembly files in .fasta format
    Alignment statistics in .csv format

Contribution

Contributions are welcome! Feel free to submit a pull request or open an issue.
License

This project is licensed under the MIT License
