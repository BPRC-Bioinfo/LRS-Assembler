# Long Read Genome Assembly with Snakemake

This Snakemake workflow assembles long reads from Oxford Nanopore (ONT) and PacBio (HiFi) datasets.
It can use individual datasets or combine both to create hybrid assemblies.
The workflow is designed for scalability and flexibility, making it ideal for genome assembly projects with varying data sources.

## Features

  * Support for ONT and/or PacBio: Assemble genomes using either sequencing platform, or combine them for a hybrid assembly.
  * Automatic Sample Detection: The pipeline automatically detects the input type and applying the appropriate assembly method.
  * Configurable Assembly Options: Users can customize provide references for alignment and comparison.
  * Reference Comparison: Align assembled sequences to reference and compute assembly statistics.

## Workflow Overview

The workflow automatically detects samples and determines the appropriate assembly method based on the available data:

  * Input Files Detection: The workflow searches for ONT .fastq.gz and PacBio .bam files in the designated directories.
  * Sample Classification: Each sample is classified into one of three categories:
        PacBio-only: Samples with only PacBio data.
        ONT-only: Samples with only ONT data.
        Hybrid (PacBio + ONT): Samples with both ONT and PacBio data.
    Assembly: Depending on the data type, the pipeline proceeds with either:
        PacBio assembly: Using tools optimized for HiFi reads.
        ONT assembly: Utilizing tools designed for long, noisy reads.
        Hybrid assembly: Combining ONT and PacBio data for higher assembly accuracy.
    Post-assembly Processing: The assembled sequences are aligned to multiple references and statistics are computed for evaluation.

Requirements

    Snakemake (>=6.0)
    Pandas for data handling
    Tools for genome assembly such as ntLink and ragtag (or any other assembly tools you wish to include).

Installation

    Clone this repository:

    git clone https://github.com/yourusername/long-read-assembly.git
    cd long-read-assembly

    Install dependencies:
        Snakemake
        Pandas
        Tools for genome assembly (e.g., ntLink, ragtag)

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

Edit the configs/run-config.yaml file to specify:

    Reference genomes: Paths to reference files used for alignment.
    Assembly parameters: Customize k-values and window sizes for assembly tools.

Here is an example run-config.yaml file:

reference:
  - path: "genomes/human/hg38.fasta"
  - path: "genomes/mouse/mm10.fasta"

Execution

Run the workflow with Snakemake:

snakemake --cores <number_of_cores>

You can also run the workflow in a cluster environment:

snakemake --cluster "qsub" --jobs 100

Rule Overview
Input Detection

Detects the input files and classifies samples based on available ONT or PacBio data:

samplesdict = {
    sample: {
        'pacbio': get_filenames(f"raws/{sample}/pacbio/*.bam"),
        'nanopore': get_filenames(f"raws/{sample}/nanopore/*.fastq.gz")
    }
    for sample in SAMPLES
}

Assembly Decision

Determines the appropriate assembly method for each sample:

def final_outcome(sample):
    pacbio_present = bool(samplesdict[sample]['pacbio'])
    nanopore_present = bool(samplesdict[sample]['nanopore'])
    if pacbio_present and nanopore_present:
        return "hybridPN"
    elif pacbio_present:
        return "pacbio"
    elif nanopore_present:
        return "nanopore"

Main Assembly Rule

The main assembly rule that processes all samples and references:

rule all:
    input:
        expand("results/{sample}/info/{sample}_{method}_denovo_ntlink_ragtag_{ref}_stats.csv",
               sample=df_final['Sample'],
               method=df_final['Assembly'],
               ref=df_final['Reference']),
    ...

Output

The output files are stored in the results/ directory. For each sample and assembly method, you will obtain:

    Assembly files in .fasta format
    Alignment statistics in .csv format

Contribution

Contributions are welcome! Feel free to submit a pull request or open an issue.
License

This project is licensed under the MIT License
