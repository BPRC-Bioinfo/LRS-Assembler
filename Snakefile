import pandas as pd

SAMPLES, NANOPORE, = glob_wildcards("raws/{sample}/nanopore/{files}.fastq.gz")
SAMPLES, PACBIO ,= glob_wildcards("raws/{sample}/pacbio/{file}.bam")

samplesdict={}
# Loop through folder to identify normal and compressed files
for sample in SAMPLES:
    samplesdict[sample]={}
    BAM, = glob_wildcards("raws/"+ sample + "/pacbio/{bam}.bam")
    samplesdict[sample]['bam'] = BAM
    FASTQ, = glob_wildcards("raws/"+ sample + "/nanopore/{fastq}.fastq.gz")
    samplesdict[sample]['fastq'] = FASTQ

# Function to extract fastq files
def extractFastq(dir,fastqType):
    return samplesdict[dir][fastqType]

rule all:
    input:
        expand("{sample}_{machine}_Mmul.sam", sample = SAMPLES, machine = ["nanopore","pacbio"])


rule pacbio_bamToFastq:
    input:
        bam = "raws/{sample}/pacbio/{pacbio}.bam",
        pbi = "raws/{sample}/pacbio/{pacbio}.bam.pbi",
    output:
        "fastq/{sample}_{pacbio}.fastq.gz"
    shell:
        """
        echo {input.bam} > {output}
        """

rule pacbio_combineFastq:
    input:
        lambda wildcards: expand("fastq/{sample}_{pacbio}.fastq.gz", sample = wildcards.sample, pacbio = extractFastq(wildcards.sample, "bam")),
    output:
        "final/{sample}_pacbio.fastq.gz"
    shell:
        """
        echo {input} > {output}
        """

rule nanopore_fastq:
    input:
        lambda wildcards: expand("raws/{sample}/nanopore/{nano}.fastq.gz", sample = wildcards.sample, nano = extractFastq(wildcards.sample, "fastq")),
    output:
        "final/{sample}_nanopore.fastq.gz"
    shell:
        """
        echo {input} > {output}
        """

rule ref_mapping:
    input:
        "final/{sample}_{machine}.fastq.gz",
    output:
        "{sample}_{machine}_Mmul.sam"
    shell:
        """
        cat {input} > {output}
        """

