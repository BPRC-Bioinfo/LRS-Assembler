# LRS-Assembly
# v 0.0.5
# By Giang Le & Jaimy


import os
import glob
import pandas as pd
from pathlib import Path

configfile: "configs/trio_t2t.yaml"

SPECIES=config['species'].replace(" ", "_")
REGIONS = config['region']

NSAMPLES = list(config["nanopore"].keys())
PSAMPLES = list(config["pacbio"].keys())

SAMPLES = list(set(NSAMPLES) & set(PSAMPLES))

print ("LRS-Assembler")
print (f"Samples detected {SAMPLES}")

REFERENCE = list(set(config["annotation_ref"].keys()))
SCAFFOLDS = list(set(config["scaffold_ref"].keys()))
print (SCAFFOLDS)

print (f"Reference detected {REFERENCE}")

wildcard_constraints:
    sample="|".join(SAMPLES),
    hap = "hap1|hap2",
    ref = "|".join(REFERENCE),
    scaffold = "|".join(SCAFFOLDS),
    region = "|".join(REGIONS)

KVALS = [24,32,40]
WVALS = [100, 250, 500]

for sample, nano_dir in config['nanopore'].items():
    print (f"For {sample} nanopore:")
    for input_dir in nano_dir:
        fastq_gz_files = glob.glob(os.path.join(input_dir, "*.fastq.gz"))
        num_files = len(fastq_gz_files)

        if num_files > 0:
            print(f"Found {num_files} fastq.gz files in {input_dir}.")
        else:
            print(f"Warning: No fastq.gz files found at {input_dir}.")

for sample, pac_dir in config['pacbio'].items():
    print (f"For {sample} pacbio:")   
    for input_dir in pac_dir:
        bam_files = glob.glob(os.path.join(input_dir, "*.bam"))
        pbi_files = glob.glob(os.path.join(input_dir, "*.bam.pbi"))
        num_bam_files = len(bam_files)
        num_pbi_files = len(pbi_files)

        if num_bam_files > 0 and num_pbi_files > 0:
            print(f"Found {num_bam_files} BAM file(s) and {num_pbi_files} PBI file(s) in {input_dir}.")
        elif num_bam_files > 0:
            print(f"Found {num_bam_files} BAM file(s), but no PBI files in {input_dir}. ")
        elif num_pbi_files > 0:
            print(f"Found {num_pbi_files} PBI file(s), but no BAM files in {input_dir}.")
        else:
            print(f"Warning: No BAM or PBI files found at {input_dir}.")

FLANK_LOCAL = {}
for region in REGIONS:

    left = REGIONS[region].get("left_flank")
    left_local = REGIONS[region].get("left_flank_local", "").strip()
    if left and left not in FLANK_LOCAL and left_local:
        FLANK_LOCAL[left] = left_local

    right = REGIONS[region].get("right_flank")
    right_local = REGIONS[region].get("right_flank_local", "").strip()
    if right and right not in FLANK_LOCAL and right_local:
        FLANK_LOCAL[right] = right_local

UNIQUE_FLANKS = set(
    [REGIONS[r]["left_flank"] for r in REGIONS] +
    [REGIONS[r]["right_flank"] for r in REGIONS]
)

def check_library(region_name, library):
    if library:
        library_path = Path(library)
        if not library_path.is_file():
            print(f"Error: Library file not found at {library} for {region_name} region.")
            sys.exit(1)
        return library
    else:
        return "no_lib"

SAMPLES = ['EAW', 'R02034']
print (SAMPLES)

rule all:
    input:
#        expand("results/{sample}/raws/{sample}_{seqMachine}_5000.fastq.gz", sample = SAMPLES, seqMachine = ['nanopore','pacbio']),
#        expand("results/{sample}/hifiasm/{sample}_hybridPN.bp.p_ctg.gfa", sample = SAMPLES),
        expand("results/{sample}/annotation/{region}/{sample}_{hap}_{species}_flanking_gene_{region}_status.csv", sample = SAMPLES, hap = ['hap1','hap2'], scaffold = SCAFFOLDS, species = SPECIES, region = REGIONS),
        expand("results/{sample}/annotation/{sample}_summary.txt", sample = SAMPLES),
        expand("results/{sample}/info/{sample}_hybridPN_scaffold_{scaffold}_infos.txt", sample = SAMPLES, scaffold = SCAFFOLDS),
        expand("results/{sample}/info/compare_{sample}_{scaffold}_figure.png", sample = SAMPLES, scaffold = SCAFFOLDS)



rule prepare_genome_reference:
    output:
        genome = "references/annotation/{ref}.fna",
    conda:
        "../envs/datasets.yaml"
    retries: 3
    params:
        genome = lambda wildcards: config["annotation_ref"][wildcards.ref]['genome'],
        ref = lambda wildcards: config["annotation_ref"][wildcards.ref]['accession_number'],
        zipfile = "{ref}_genome.zip",
    shell:
        """
        if [ -z "{params.genome}" ] || [ ! -f "{params.genome}" ]; then
            echo "Genome file not found locally, downloading from NCBI"
            datasets download genome accession {params.ref} --filename {params.zipfile} --include genome || (echo "Warning: download genome failed. Check accession or prepare the genome locally.")
            unzip {params.zipfile} -d "{params.ref}_fna"
            mv {params.ref}_fna/ncbi_dataset/data/{params.ref}/{params.ref}*.fna {output.genome}
            rm -r {params.zipfile} {params.ref}_fna
        else
            echo "Using local genome file {params.genome}"
            cp {params.genome} {output.genome}
        fi
        echo "Complete prepare_genome_reference"
        """

rule prepare_gff_reference:
    output:
        gff = "references/annotation/{ref}.gff",
    conda:
        "../envs/datasets.yaml"
    retries: 3
    params:
        ref = lambda wildcards: config["annotation_ref"][wildcards.ref]['accession_number'],
        gff = lambda wildcards: config["annotation_ref"][wildcards.ref]['gff'],
        zipfile = "{ref}_gff.zip",
    shell:
        """
        if [ -z "{params.gff}" ] || [ ! -f "{params.gff}" ]; then
            echo "Gff file not found locally, downloading from NCBI"
            datasets download genome accession {params.ref} --filename {params.zipfile} --include gff3 || (echo "Warning: download gff failed. Check accession or prepare locally.")
            unzip {params.zipfile} -d "{params.ref}_gff"
            mv {params.ref}_gff/ncbi_dataset/data/{params.ref}/genomic.gff {output.gff}
            rm -r {params.zipfile} {params.ref}_gff
        else
            echo "Using local genome file {params.gff}"
            cp {params.gff} {output.gff}
        fi
        echo "Complete prepare_gff_reference"
        """

rule prepare_scaffold_reference:
    output:
        genome = "references/scaffolds/{scaffold}.fna",
    conda:
        "../envs/datasets.yaml"
    retries: 3
    params:
        genome = lambda wildcards: config["scaffold_ref"][wildcards.scaffold]['genome'],
        ref = lambda wildcards: config["scaffold_ref"][wildcards.scaffold]['accession_number'],
        zipfile = "{scaffold}_genome.zip",
    shell:
        """
        if [ -z "{params.genome}" ] || [ ! -f "{params.genome}" ]; then
            echo "Genome file not found locally, downloading from NCBI"
            datasets download genome accession {params.ref} --filename {params.zipfile} --include genome || (echo "Warning: download genome failed. Check accession or prepare the genome locally.")
            unzip {params.zipfile} -d "{params.ref}_fna"
            mv {params.ref}_fna/ncbi_dataset/data/{params.ref}/{params.ref}*.fna {output.genome}
            rm -r {params.zipfile} {params.ref}_fna
        else
            echo "Using local genome file {params.genome}"
            cp {params.genome} {output.genome}
            echo "Prepare the info file locally for your custom reference genome."
        fi
        echo "Complete prepare_genome_reference"
        """

rule prepare_chromosome_reference:
    output:
        info = "references/scaffolds/{scaffold}.info",
    retries: 3
    conda:
        "../envs/datasets.yaml"
    params:
        ref = lambda wildcards: config["scaffold_ref"][wildcards.scaffold]['accession_number'],
        info = lambda wildcards: config["scaffold_ref"][wildcards.scaffold]['chr_info'],
        zipfile = "{scaffold}_info.zip",
    shell:
        """
        if [ -z "{params.info}" ] || [ ! -f "{params.info}" ]; then
            echo "Info file not found locally, downloading from NCBI"
            datasets download genome accession {params.ref} --filename {params.zipfile} --include seq-report || (echo "Warning: download info failed. Check accession or prepare locally.")
            unzip {params.zipfile} -d "{params.ref}_info"
            python scripts/accession_to_chromosome.py {params.ref}_info/ncbi_dataset/data/{params.ref}/sequence_report.jsonl | awk '{{print $1"\t"$2}}' > {output.info}
            rm -r {params.zipfile} {params.ref}_info
        else
            echo "Using local genome file {params.info}"
            cp {params.info} {output.info}
        fi
        echo "Complete prepare_chromosome_reference"
        """

rule prepare_chromosome1_reference:
    output:
        info = "references/{ref}.info",
    retries: 3
    conda:
        "../envs/datasets.yaml"
    params:
        ref = lambda wildcards: config["reference"][wildcards.ref]['accession_number'],
        info = lambda wildcards: config["reference"][wildcards.ref]['chr_info'],
        zipfile = "{ref}_info.zip",
    shell:
        """
        if [ -z "{params.info}" ] || [ ! -f "{params.info}" ]; then
            echo "Info file not found locally, downloading from NCBI"
            datasets download genome accession {params.ref} --filename {params.zipfile} --include seq-report || (echo "Warning: download info failed. Check accession or prepare locally.")
            unzip {params.zipfile} -d "{params.ref}_info"
            python scripts/json_parse.py {params.ref}_info/ncbi_dataset/data/{params.ref}/sequence_report.jsonl | sed '/accession/d;/MT/d' | awk '{{print $1",@_chr"$2}}' > {output.info}
            rm -r {params.zipfile} {params.ref}_info
        else
            echo "Using local genome file {params.info}"
            cp {params.info} {output.info}
        fi
        echo "Complete prepare_chromosome_reference"
        """


## Process raw files
rule combine_filter_nanopore_fastq:
    input:
        nanopore_files = lambda wildcards: [f for dir in config["nanopore"][wildcards.sample]
                                            for f in glob.glob(f"{dir}/*.fastq.gz")]
    output:
        "results/{sample}/raws/{sample}_nanopore_5000.fastq.gz"
    threads: 5
    benchmark:
        "results/{sample}/benchmarks/01_{sample}_nanopore_combine_filter_fastq.bench"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        if [ -z "{input}" ]; then
            echo "No nanopore fastq files found"
            exit
        fi
        for i in {input}; do
            cat $i
        done | seqkit sana | seqkit seq -m 5000 -Q 9 | seqkit rmdup -s -o {output}
        echo "The rule combine_filter_nanopore_fastq complete for sample {wildcards.sample}"
        """

rule combine_pacbio_fastq:
    input:
        pacbio_bam = lambda wildcards: [f for dir in config["pacbio"][wildcards.sample]
                                            for f in glob.glob(f"{dir}/*.bam")],
        pacbio_pbi = lambda wildcards: [f for dir in config["pacbio"][wildcards.sample]
                                            for f in glob.glob(f"{dir}/*.bam.pbi")]
    output:
        temp("fastq/{sample}_pacbio_raw.fastq.gz")
    threads: 8
    benchmark:
        "results/{sample}/benchmarks/01_{sample}_pacbio_combine_fastq.bench"
    conda:
        "../envs/bam2fastx.yaml"
    shell:
        """
        if [ -z "{input.pacbio_bam}" ]; then
            echo "No pacbio bam files found"
            exit
        fi
        bam2fastq -o fastq/{wildcards.sample}_pacbio_raw {input.pacbio_bam}
        echo "The rule combine_pacbio_fastq complete for sample {wildcards.sample}"
        """

rule filter_pacbio_fastq:
    input:
        rules.combine_pacbio_fastq.output
    output:
        "results/{sample}/raws/{sample}_pacbio_5000.fastq.gz"
    threads: 5
    benchmark:
        "results/{sample}/benchmarks/01_{sample}_pacbio_filter_fastq.bench"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit sana {input} | seqkit seq -m 5000 | seqkit rmdup -s -o {output}
        echo "The rule filter_pacbio_fastq complete for sample {wildcards.sample}"
        """

rule filtered_reads_stats:
    input:
        filtered = ancient("results/{sample}/raws/{sample}_{seqMachine}_5000.fastq.gz")
    output:
        filtered = "results/{sample}/info/{sample}_{seqMachine}_stat_filtered.txt",
    wildcard_constraints:
        seqMachine="(nanopore|pacbio)"
    threads: 5
    conda:
        "../envs/seqkit.yaml"
    benchmark:
        "results/{sample}/benchmarks/01_{sample}_{seqMachine}_fastq_filter_stats.bench"
    threads: 5
    shell:
        """
        seqkit stats --threads {threads} {input.filtered} > {output.filtered}
        echo "The rule filter_read_stats complete for sample {wildcards.sample} {wildcards.seqMachine}"
        """

## Genome assembling
rule hifiasmUL_de_novo:
    input:
        nano = "results/{sample}/raws/{sample}_nanopore_5000.fastq.gz",
        pac = "results/{sample}/raws/{sample}_pacbio_5000.fastq.gz",
    output:
        prima = "results/{sample}/hifiasm/{sample}_hybridPN.bp.p_ctg.gfa",
        hap1 = "results/{sample}/hifiasm/{sample}_hybridPN.bp.hap1.p_ctg.gfa",
        hap2 = "results/{sample}/hifiasm/{sample}_hybridPN.bp.hap2.p_ctg.gfa",
    conda:
        "../envs/hifiasm.yaml"
    log:
        "results/{sample}/logs/01_{sample}_hybridPN_hifiasm_denovo.log"
    threads: 20
    params:
        outdir = "results/{sample}/hifiasm/{sample}_hybridPN",
        hifiasm_settings = config['hifiasm_settings']
    benchmark:
        "results/{sample}/benchmarks/02_{sample}_hybridPN_hifiasm_denovo_assembly.bench"
    shell:
        """
        if [[ -z "{params.hifiasm_settings}" ]]; then
            hifiasm -o {params.outdir} -t {threads} --ul {input.nano} {input.pac} 2> {log}
        else
            hifiasm -o {params.outdir} -t {threads} --ul {input.nano} {params.hifiasm_settings} {input.pac} 2> {log}
        fi
        echo "The rule hifiasmUL_de_novo complete for sample {wildcards.sample}"
        """

rule gfa_to_fasta_primary:
    input:
        prima = ancient(rules.hifiasmUL_de_novo.output['prima']),
    output:
        prima = "results/{sample}/midway/{sample}_hybridPN_denovo.fa",
    benchmark:
        "results/{sample}/benchmarks/02_{sample}_primary_assembly_gfa_to_fasta.bench"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input.prima} | sed "s/ptg/{wildcards.sample}_/g" > {output.prima}
        """

rule primary_contigs_quast:
    input:
        rules.gfa_to_fasta_primary.output.prima,
    output:
        "results/{sample}/info/{sample}_primary_contig_report.tsv"
    conda:
        "../envs/quast.yaml"
    benchmark:
        "results/{sample}/benchmarks/02_{sample}_hybridPN_quast.bench"
    params:
        "results/{sample}/info/{sample}_hybridPN_denovo_quast/"
    shell:
        """
        quast.py {input} -o {params}
        mv {params}"report.tsv" {output}
        rm -r {params}
        echo "The rule contig_quast complete for sample {wildcards.sample}"
        """

rule gfa_to_fasta_haps:
    input:
        hap1 = ancient(rules.hifiasmUL_de_novo.output['hap1']),
        hap2 = ancient(rules.hifiasmUL_de_novo.output['hap2']),
    output:
        hap1 = "results/{sample}/midway/{sample}_hybridPN_denovo_hap1.fa",
        hap2 = "results/{sample}/midway/{sample}_hybridPN_denovo_hap2.fa",
    benchmark:
        "results/{sample}/benchmarks/02_{sample}_haplotype_assemblies_gfa_to_fasta.bench"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input.hap1} | sed "s/h1tg/{wildcards.sample}_hap1_/g" > {output.hap1}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap2} | sed "s/h2tg/{wildcards.sample}_hap2_/g" > {output.hap2}
        """

rule primary_genome_scaffold:
    input:
        ref = ancient(rules.prepare_scaffold_reference.output),
        fa = rules.gfa_to_fasta_primary.output.prima,
    output:
        temp("results/{sample}/scaffolds/{sample}_hybridPN_scaffold_{scaffold}.fa")
    params:
        "results/{sample}/midway/ragtag/{scaffold}/scaffold/"
    threads: 10
    benchmark:
        "results/{sample}/benchmarks/04_{sample}_hybridPN_{scaffold}_ragtag-scaffold.bench"
    conda:
        "../envs/ragtag.yaml"
    shell:
        """
        if [[ {wildcards.scaffold} != "none" ]]; then
            mkdir -p {params}
            ragtag.py scaffold -t {threads} -o {params} -u -r -g 1 -m 2000000000 {input.ref} {input.fa}
            sed 's/_RagTag//g' {params}"ragtag.scaffold.fasta" > {output}
            rm -r {params}
        else
            cp {input.fa} {output}
        fi
        """

rule rename_primary_scaffold:
    input:
        acc = rules.prepare_chromosome_reference.output,
        fa = rules.primary_genome_scaffold.output,
    output:
        fa = "results/{sample}/scaffolds/{sample}_hybridPN_scaffold_{scaffold}.fasta",
    conda:
        "../envs/rename.yaml"
    shell:
        """
        python scripts/rename_accession.py --ac_file {input.acc} --input_fasta {input.fa} --output_fasta {output} --header_pattern "{wildcards.sample}" 
        """

rule primary_scaffold_quast:
    input:
        rules.rename_primary_scaffold.output,
    output:
        "results/{sample}/info/{sample}_scaffold_{scaffold}_quast_report.tsv"
    conda:
        "../envs/quast.yaml"
    benchmark:
        "results/{sample}/benchmarks/06_{sample}_hybridPN_{scaffold}_quast_final.bench"
    params:
        "results/{sample}/scaffolds/{sample}_scaffold_{scaffold}_quast/"
    shell:
        """
        quast.py {input} -o {params}
        mv {params}"report.tsv" {output}
        rm -r {params}
        echo "Complete quast analysis of sample {wildcards.sample} scaffolds"
        """

checkpoint primary_scaffold_busco:
    input:
        rules.rename_primary_scaffold.output,
    output:
        busco_dir = directory("results/{sample}/info/{sample}_scaffold_{scaffold}_busco/"),
    conda:
        "../envs/busco.yaml"
    log:
        "results/{sample}/logs/03_{sample}_scaffold_{scaffold}_busco.log"
    benchmark:
        "results/{sample}/benchmarks/06_{sample}_hybridPN_{scaffold}_busco.bench"
    params:
        busco = config["busco"],       
    threads: 20
    shell:
        """
        busco -m genome -i {input} -o {output.busco_dir} -l {params} -c {threads} -f 2> {log}
        echo "Complete BUSCO analysis of sample {wildcards.sample} scaffolds"
        """

def primary_scaffold_busco_output(wildcards):
    checkpoint_output = checkpoints.primary_scaffold_busco.get(**wildcards).output[0]
    return expand("results/{sample}/info/{sample}_scaffold_{scaffold}_busco/short_summary.{i}_busco.txt",
                sample = wildcards.sample,
                scaffold=wildcards.scaffold,
                i=glob_wildcards(os.path.join(checkpoint_output, "short_summary.{i}_busco.txt")).i)

checkpoint reference_busco:
    input:
        rules.prepare_scaffold_reference.output,
    output:
        directory("references/scaffolds/{scaffold}_busco/"),
    conda:
        "../envs/busco.yaml"
    params:
        busco = config["busco"],
    threads: 20
    shell:
        """
        busco -m genome -i {input} -o {output} -l {params.busco} -c {threads} -f 
        """

def reference_busco_output(wildcards):
    checkpoint_output = checkpoints.reference_busco.get(**wildcards).output[0]
    return expand("references/scaffolds/{scaffold}_busco/short_summary.{i}_busco.txt",
                scaffold=wildcards.scaffold,
                i=glob_wildcards(os.path.join(checkpoint_output, "short_summary.{i}_busco.txt")).i)

rule busco_directory:
    input:
        ref = reference_busco_output,
        scaffold = primary_scaffold_busco_output
    output:
        directory("results/{sample}/info/compare_{sample}_{scaffold}_busco/")
    shell:
        """
        mkdir -p {output}
        cp {input} {output}
        """

rule busco_final_figure:
    input:
        file = "results/{sample}/info/{sample}_vs_{scaffold}_busco.txt",
        outdir ="results/{sample}/info/compare_{sample}_{scaffold}_busco/"
    output:
        "results/{sample}/info/compare_{sample}_{scaffold}_figure.png"
    conda:
        "../envs/busco.yaml"
    benchmark:
        "results/{sample}/benchmarks/06_{sample}_hybridPN_{scaffold}_busco_image.bench"
    shell:
        """
        generate_plot.py --working_directory {input.outdir}
        mv {input.outdir}/busco_figure.png {output}
        """

rule busco_compare:
    input:
        ref = reference_busco_output,
        scaffold = primary_scaffold_busco_output
    output:
        "results/{sample}/info/{sample}_vs_{scaffold}_busco.txt"
    shell:
        """
        custom_array=("Complete" "Complete_and_single-copy" "Complete_and_duplicated" "Fragmented" "Missing" "Total")
        (echo "Description, {wildcards.scaffold}, {wildcards.sample}"; paste <(grep "Result" -A 8 {input.ref} | tail -6 | awk '{{print $1}}') <(grep "Result" -A 8 {input.scaffold} | tail -6 | awk '{{print $1}}')  | awk -v arr="${{custom_array[*]}}" 'BEGIN {{split(arr, a, " ")}} {{printf "%s, %s, %s\\n", a[NR], $1, $2}}') > {output}
        """

rule primary_chromosomes_only:
    input:
        rules.rename_primary_scaffold.output
    output:
        fa = "results/{sample}/scaffolds/{sample}_hybridPN_scaffold_{scaffold}_chrs.fasta",
        names = temp("results/{sample}/scaffolds/{sample}_hybridPN_scaffold_{scaffold}_chrs.names"),
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        grep "chr" {input} | grep -v "MT" | sed 's/>//g' > {output.names}
        seqkit grep -f {output.names} {input} > {output.fa}
        """

rule primary_chromosomes_Ngaps:
    input:
        fa = rules.primary_chromosomes_only.output.fa,
    output:
        "results/{sample}/info/{sample}_hybridPN_scaffold_{scaffold}_chr_gaps.txt",
    shell:
        """
        python scripts/N-detect.py {input.fa} > {output}
        """

rule primary_chromosomes_unplaced:
    input:
        fa = rules.rename_primary_scaffold.output.fa,
        names = rules.primary_chromosomes_only.output.names,
    output:
        "results/{sample}/info/{sample}_hybridPN_scaffold_{scaffold}_unplaced.txt",
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit grep -v -f {input.names} {input.fa} | seqkit stats > {output}
        """


rule assembly_infos:
    input:
        contig = rules.primary_contigs_quast.output,
        unplaced = rules.primary_chromosomes_unplaced.output,
        Ngaps = rules.primary_chromosomes_Ngaps.output,
        scaffold = rules.primary_scaffold_quast.output,
        busco = primary_scaffold_busco_output
    output:
        "results/{sample}/info/{sample}_hybridPN_scaffold_{scaffold}_infos.txt"
    shell:
        """
        # Genome size
        grep "Total length" {input.contig}| grep " 0 bp" | awk -F'\t' '{{print $2}}' > {output}
        # Unplaced length
        grep "-" {input.unplaced} | awk '{{print $5}}' | tr -d ',' >> {output}
        # Contigs
        grep "contigs" {input.contig} | grep " 0 bp"  | awk -F'\t' '{{print $2}}' >> {output}
        # N50 contig stats
        grep "N50" {input.contig} | awk -F'\t' '{{print $2}}' >> {output}
        # Gaps
        wc -l {input.Ngaps} | awk '{{print $1}}' >> {output}
        # Gap length
        sed 's/.*://g' {input.Ngaps} | awk 'SUM+=$0; END {{print SUM}}' | tail -1 >> {output}
        # Scaffolds
        grep "contigs" {input.scaffold} | grep " 0 bp"  | awk -F'\t' '{{print $2}}' >> {output}
        # N50 scaffold stats
        grep "N50" {input.scaffold} | awk -F'\t' '{{print $2}}' >> {output}
        # Busco
        total_busco=$(grep "Total BUSCO" {input.busco} | awk '{{print $1}}')
        grep "(C)" {input.busco} | awk -v total="$total_busco" '{{print $1/total * 100}}' >> {output}
        grep "(C)" {input.busco} | awk '{{print $1}}' >> {output}
        """


## Fixing up to here

rule haplotype_scaffold:
    input:
        ref = rules.prepare_scaffold_reference.output,
        fa = "results/{sample}/midway/{sample}_hybridPN_denovo_{hap}.fa"
    output:
        temp("results/{sample}/midway/scaffolds/{sample}_hybridPN_{hap}_{scaffold}.fa")
    params:
        "results/{sample}/midway/scaffolds/{scaffold}/{hap}/"
    threads: 10
    benchmark:
        "results/{sample}/benchmarks/04_{sample}_hybridPN_{hap}_{scaffold}.bench"
    conda:
        "../envs/ragtag.yaml"
    shell:
        """
        if [[ {wildcards.scaffold} != "none" ]]; then
            mkdir -p {params}
            ragtag.py scaffold -t {threads} -o {params} -u -r -g 1 -m 2000000000 {input.ref} {input.fa}
            sed 's/_RagTag//g' {params}"/ragtag.scaffold.fasta" > {output}
            rm -r {params}
        else
            cp {input.fa} {output}
        fi
        """

rule rename_haplo_scaffolds:
    input:
        acc = rules.prepare_chromosome_reference.output,
        fa = rules.haplotype_scaffold.output,
    output:
        fa = temp("results/{sample}/midway/scaffolds/{sample}_hybridPN_{hap}_{scaffold}.fasta"),
    conda:
        "../envs/rename.yaml"
    benchmark:
        "results/{sample}/benchmarks/04_{sample}_hibridPN_{scaffold}_{hap}_rename_scaffolds.bench"
    shell:
        """
        python scripts/rename_accession.py --ac_file {input.acc} --input_fasta {input.fa} --output_fasta {output} --header_pattern "{wildcards.sample}_{wildcards.hap}" 
        """

## Genome annotation
rule scaffolded_haplotypes:
    input:
        expand("results/{{sample}}/midway/scaffolds/{{sample}}_hybridPN_{{hap}}_{scaffold}.fasta", scaffold = SCAFFOLDS),
    output:
        "results/{sample}/scaffolds/{sample}_hybridPN_{hap}_scaffolded.fa"
    shell:
        """
        cat {input} > {output}
        """

## Working up to here


## Annotation analysis
rule prepare_flanking_genes:
    output:
        "flanking_genes/{species}/{flank}.fasta"
    params:
        local_file=lambda wc: FLANK_LOCAL.get(wc.flank, ""),
        species=config['species']
    retries: 3
    conda:
        "../envs/flanking.yaml"
    shell:
        r"""
        if [ -n "{params.local_file}" ] && [ -f "{params.local_file}" ]; then
            cp {params.local_file} {output}
        else
            MINWAIT=3
            MAXWAIT=6
            sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))
            datasets download gene symbol {wildcards.flank} --taxon "{params.species}" --include gene --filename {wildcards.flank}.zip
            unzip {wildcards.flank}.zip -d {wildcards.species}_{wildcards.flank}
            cat {wildcards.species}_{wildcards.flank}/ncbi_dataset/data/gene.fna > {output}
            rm -r {wildcards.species}_{wildcards.flank} {wildcards.flank}.zip
        fi
        """

rule flanking_genes_in_region:
    input:
        flanks = lambda wc: [
            f"flanking_genes/{SPECIES}/{REGIONS[wc.region][key].strip()}.fasta"
            for key in ["left_flank", "right_flank"]
            if REGIONS[wc.region].get(key, "").strip() 
        ]
    output:
        fasta = "regions/{species}/{region}.fasta",
    shell:
        """
        cat {input.flanks} | sed -E "s/ \\[.*//g;s/ /_/g" > {output}
        """

rule region_or_telomere:
    input:
        rules.flanking_genes_in_region.output
    output:
        "regions/{species}/{region}.info"
    params:
        left = lambda wc: REGIONS[wc.region].get("left_flank", "").strip(),
        right = lambda wc: REGIONS[wc.region].get("right_flank", "").strip()
    shell:
        """
        if [ -n "{params.left}" ]; then
            echo "left: {params.left}" >> {output}
        fi
        if [ -n "{params.right}" ]; then
            echo "right: {params.right}" >> {output}
        fi
        """

rule mapping_flanking_genes:
    input:
        fa = "results/{sample}/scaffolds/{sample}_hybridPN_{hap}_scaffolded.fa",
        flanks = "regions/{species}/{region}.fasta",
    output:
        sam = "map_temp/{species}/{sample}/{region}/{sample}_{hap}_flank_genes_{region}.sam",
        flank_region = "map_temp/{species}/{sample}/{region}/{sample}_{hap}_flank_genes_{region}.txt"
    log:
        "map_temp/{species}/{sample}/{region}/logs/{sample}_{hap}_{region}_flanking.log"
    benchmark:
        "results/{sample}/benchmarks/07_{sample}_{hap}_{species}_{region}_1_map_flanking_genes.bench"
    threads: 6
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax asm5 -t {threads} {input.fa} {input.flanks} > {output.sam} 2> {log}
        cut -f1-5 {output.sam} | grep -v "@" | awk '$3!="*"' > {output.flank_region}
        """

checkpoint check_flanking_genes:
    input:
        loc = rules.mapping_flanking_genes.output['flank_region'],
        info = rules.region_or_telomere.output,
    output:
        "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_flanking_gene_{region}_status.csv",
    shell:
        """
        python scripts/check_flanking_genes.py -i {input.info} -l {input.loc} -o {output}
        """

rule region_sam_to_bed:
    input:
        sam = rules.mapping_flanking_genes.output['sam'],
    output:
        bam = temp("map_temp/{species}/{sample}/{region}/{sample}_{hap}_{species}_flank_genes_{region}.bam"),
        sort = temp("map_temp/{species}/{sample}/{region}/{sample}_{hap}_{species}_flank_genes_{region}_sort.bam"),
        bed = "results/{sample}/annotation/{region}/flanking/{sample}_{hap}_{species}_flanking_genes_{region}.bed",
    threads: 2
    benchmark:
        "results/{sample}/benchmarks/07_{sample}_{hap}_{species}_{region}_3_sam_bed.bench"
    conda:
        "../envs/sambed.yaml"
    shell:
        """
        samtools view -@ {threads} -bh {input.sam} > {output.bam}
        samtools sort -@ {threads} -o {output.sort} {output.bam}
        bedtools bamtobed -i {output.sort} > {output.bed}
        """

rule intact_region_coordinates:
    input:
        bed = rules.region_sam_to_bed.output['bed'],
        status = rules.check_flanking_genes.output,
    output:
        "results/{sample}/annotation/{region}/flanking/{sample}_{hap}_{species}_{region}_intact.bed",
    shell:
        """
        python scripts/bam2bed.py -s {input.status} -b {input.bed} -r {wildcards.region} -o {output}
        """

rule intact_region_contig_extraction:
    input:
        bed = rules.intact_region_coordinates.output,
        fa = rules.scaffolded_haplotypes.output,
    output:
        "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}_intacted.fa",
    params:
        "results/{sample}/scaffolds/{sample}_hybridPN_{hap}_scaffolded.fa.fai"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools getfasta -fi {input.fa} -bed {input.bed} -s | sed 's/:.*/_{wildcards.region}/g' > {output}
        if  [[ -f {params} ]]; then
            rm {params}
        fi
        """

checkpoint libraries_format:
    output:
        directory("map_temp/{species}/libraries/{region}/"),
    params:
        cDNA_in = lambda wc: check_library(
            wc.region,
            config['region'][wc.region].get('cDNA_library', False)
        ),
        gDNA_in = lambda wc: check_library(
            wc.region,
            config['region'][wc.region].get('gDNA_library', False)
        ),
        protein_in = lambda wc: check_library(
            wc.region,
            config['region'][wc.region].get('protein_library', False)
        ),
        cDNA_out = "map_temp/{species}/libraries/{region}/cDNA.fasta",
        gDNA_out = "map_temp/{species}/libraries/{region}/gDNA.fasta",
        protein_out = "map_temp/{species}/libraries/{region}/protein.fasta",
    conda:
        "../envs/rename.yaml"
    shell:
        """
        mkdir -p {output}
        if [ {params.cDNA_in} != "no_lib" ] ; then
            python scripts/lib_format.py {params.cDNA_in} | seqkit rmdup -s | sed 's/\\r//g' > {params.cDNA_out}
        fi

        if [ {params.gDNA_in} != "no_lib" ]; then
            python scripts/lib_format.py {params.gDNA_in} | seqkit rmdup -s | sed 's/\\r//g' > {params.gDNA_out}
        fi

        if [ {params.protein_in} != "no_lib" ]; then
            python scripts/lib_format.py {params.protein_in} | seqkit rmdup -s | sed 's/\\r//g' > {params.protein_out}
        fi
        """

rule roi_cDNA_mapping:
    input:
        fa = "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}_{roi}.fa",
        library = "map_temp/{species}/libraries/{region}/cDNA.fasta"
    output:
        "map_temp/{species}/{sample}/{region}/analysis/{sample}_{hap}_{species}_{region}_{roi}_cDNA.psl",
    threads: 10
    params:
        "map_temp/{species}/{sample}/{region}/analysis/cDNA/"
    conda:
        "../envs/gmap.yaml"
    shell:
        """
        gmap_build -D {params} -d {wildcards.sample}_{wildcards.hap}_{wildcards.region} {input.fa}
        gmap -t {threads} -D {params} -d {wildcards.sample}_{wildcards.hap}_{wildcards.region} {input.library} > {output}
        """

rule roi_parcing_gmap:
    input:
        rules.roi_cDNA_mapping.output,
    output:
        "results/{sample}/annotation/{region}/candidates/{sample}_{hap}_{species}_{region}_{roi}_candidates_cDNA.csv"
    log:
        "map_temp/{species}/{sample}/{region}/logs/{sample}_{hap}_{species}_{region}_{roi}_cDNA.log"
    shell:
        """
        python scripts/gmap_parse.py {input} {output}
        """

rule roi_gDNA_mapping:
    input:
        fa = "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}_{roi}.fa",
        library = "map_temp/{species}/libraries/{region}/gDNA.fasta"
    output:
        "map_temp/{species}/{sample}/{region}/analysis/{sample}_{hap}_{species}_{region}_{roi}_gDNA.paf",
    threads: 6
    log:
        "map_temp/{species}/{sample}/{region}/logs/{sample}_{hap}_{species}_{region}_{roi}_gDNA.log"
    benchmark:
        "results/{sample}/benchmarks/07_{sample}_{hap}_{species}_{species}_{region}_6_lib_maps_{roi}.bench"
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -x asm5 -c -t {threads} {input.fa} {input.library} > {output} 2> {log}
        """

rule roi_cigar_info:
    input:
        rules.roi_gDNA_mapping.output,
    output:
        "results/{sample}/annotation/{region}/candidates/{sample}_{hap}_{species}_{region}_{roi}_coords_paf_gDNA.txt"
    benchmark:
        "results/{sample}/benchmarks/07_{sample}_{hap}_{species}_{region}_{roi}_7_cigar_process.bench"
    shell:
        """
        python scripts/cigar_digest.py {input} {output}
        """

rule roi_gene_blast:
    input:
        paf = "results/{sample}/annotation/{region}/candidates/{sample}_{hap}_{species}_{region}_{roi}_coords_paf_gDNA.txt",
        fa = "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}_{roi}.fa",
        lib = "map_temp/{species}/libraries/{region}/gDNA.fasta",
    output:
        blast = "map_temp/{species}/{sample}/{region}/analysis/{sample}_{hap}_{species}_{region}_{roi}_gDNA.blast",
        fai = temp("results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}_{roi}.fa.fai"),
    log:
        "results/{sample}/annotation/{region}/logs/{sample}_{hap}_{species}_{region}_{roi}_blast_gDNA.log"
    benchmark:
        "results/{sample}/benchmarks/07_{sample}_{hap}_{species}_{region}_8_gene_blast_{roi}.bench"
    conda:
        "../envs/bbstools.yaml"
    threads: 5
    shell:
        """
        python scripts/paf_blast.py -p {input.paf} -r {input.lib} -f {input.fa} -o {output.blast} &> {log}
        """

rule roi_join_blast_paf:
    input:
        info = rules.roi_cigar_info.output,
        blast = rules.roi_gene_blast.output.blast,
    output:
        "results/{sample}/annotation/{region}/candidates/{sample}_{hap}_{species}_{region}_{roi}_candidates_gDNA.csv"
    shell:
        """
        python scripts/blast_extraction.py -i {input.info} -b {input.blast} -o {output}
        """

rule roi_protein_mapping:
    input:
        fa = "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}_{roi}.fa",
        library = "map_temp/{species}/libraries/{region}/protein.fasta"
    output:
        "map_temp/{species}/{sample}/{region}/analysis/{sample}_{hap}_{species}_{region}_{roi}_protein.paf",
    threads: 6
    log:
        "map_temp/{species}/{sample}/{region}/logs/{sample}_{hap}_{species}_{region}_{roi}_protein.log"
    conda:
        "../envs/miniprot.yaml"
    shell:
        """
        miniprot -t {threads} {input.fa} {input.library} --aln --trans > {output} 2> {log}
        """

rule roi_protein_extraction:
    input:
        rules.roi_protein_mapping.output
    output:
        "results/{sample}/annotation/{region}/candidates/{sample}_{hap}_{species}_{region}_{roi}_candidates_protein.csv"
    shell:
        """
        python scripts/miniprot_parse2.py {input} {output}
        """

rule roi_annotation_filter:
    input:
        "results/{sample}/annotation/{region}/candidates/{sample}_{hap}_{species}_{region}_{roi}_candidates_{i}.csv"
    output:
        multiext("results/{sample}/annotation/{region}/single_libs/{sample}_{hap}_{species}_{region}_{roi}_final_{i}", ".csv", ".bed")
    params:
        "results/{sample}/annotation/{region}/single_libs/{sample}_{hap}_{species}_{region}_{roi}_final_{i}"
    shell:
        """
        python scripts/filter_annotation2.py -i {input} -l {wildcards.i} -o {params}
        """

## Annotation missing
rule missing_raws_to_reference:
    input:
        fq = "results/{sample}/raws/{sample}_{seqMachine}_5000.fastq.gz", 
        ref = "references/scaffolds/{scaffold}.fna",
    output:
        sam = temp("map_temp/{species}/{sample}/{scaffold}_{seqMachine}.sam"),
    threads: 6
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        if [[ {wildcards.seqMachine} == "nanopore" ]]; then
            minimap2 -ax map-ont -t {threads} {input.ref} {input.fq} > {output.sam} 
        else
            minimap2 -ax map-hifi -t {threads} {input.ref} {input.fq} > {output.sam}
        fi
        """

rule missing_raws_to_reference_sam_to_bam:
    input:
        rules.missing_raws_to_reference.output,
    output:
        bam = temp("map_temp/{species}/{sample}/{scaffold}_{seqMachine}.bam"),
    threads: 6
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -bh {input} > {output.bam}
        """

rule missing_raws_to_reference_sort:
    input:
        rules.missing_raws_to_reference_sam_to_bam.output,
    output:
        sorted = temp("map_temp/{species}/{sample}/{scaffold}_{seqMachine}_sorted.bam"),
        sorted_bai = temp("map_temp/{species}/{sample}/{scaffold}_{seqMachine}_sorted.bam.bai"),
    threads: 8
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -o {output.sorted} {input}
        samtools index {output.sorted}
        """

rule mapping_roi_flanking_genes_to_reference_genome:
    input:
        fa = rules.prepare_scaffold_reference.output, 
        flanks = rules.flanking_genes_in_region.output,
    output:
        sam = temp("map_temp/{species}/{sample}/{region}/reference_{scaffold}_{region}.sam"),
        roi_region = "map_temp/{species}/{sample}/{region}/reference_{scaffold}_{region}.txt"
    threads: 6
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax asm5 -t {threads} {input.fa} {input.flanks} > {output.sam}
        cut -f1-5 {output.sam} | grep -v "@" | awk '$3!="*"' > {output.roi_region}
        """

rule missing_reference_chromosomes_extract:
    input:
        contig = rules.mapping_roi_flanking_genes_to_reference_genome.output.roi_region,
        fa = rules.prepare_scaffold_reference.output,
    output:
        temp = temp("map_temp/{species}/{sample}/{region}/reference_{scaffold}_{region}.chrs"),
        fa = "map_temp/{species}/{sample}/{region}/reference_{scaffold}_{region}.fasta"
    threads: 6
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        awk '{{print $3}}' {input.contig} | awk '!s[$0]++' > {output.temp}
        seqkit grep -f {output.temp} {input.fa} > {output.fa}
        """

rule identify_roi_scaffolds_from_reference_genome:
    input:
        region = rules.mapping_roi_flanking_genes_to_reference_genome.output.roi_region,
        raws_bam = rules.missing_raws_to_reference_sort.output.sorted,
        sorted_bai = rules.missing_raws_to_reference_sort.output.sorted_bai
    output:
        temp("map_temp/{species}/{sample}/{region}/{sample}_{hap}_{scaffold}_{region}_{seqMachine}.bam")
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        # What if on multiple scaffolds?
        awk '!s[$3]++ {{print $3}}' {input.region} | for i in `cat`; do
            samtools view -b {input.raws_bam} $i > {output}
        done
        """

rule missing_roi_raws:
    input:
        rules.identify_roi_scaffolds_from_reference_genome.output,
    output:
        "map_temp/{species}/{sample}/{region}/{sample}_{hap}_{scaffold}_{region}_{seqMachine}.fastq.gz"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools fastq -0 {output} {input}
        """

rule hifiasmUL_reference_roi:
    input:
        nano = "map_temp/{species}/{sample}/{region}/{sample}_{hap}_{scaffold}_{region}_nanopore.fastq.gz",
        pac = "map_temp/{species}/{sample}/{region}/{sample}_{hap}_{scaffold}_{region}_pacbio.fastq.gz",
    output:
        prima = "map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_{scaffold}_roi.bp.p_ctg.gfa",
        hap1 = "map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_{scaffold}_roi.bp.hap1.p_ctg.gfa",
        hap2 = "map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_{scaffold}_roi.bp.hap2.p_ctg.gfa",
    conda:
        "../envs/hifiasm.yaml"
    threads: 20
    params:
        "map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_{scaffold}_roi"
    shell:
        """
        hifiasm -o {params} -t {threads} --ul {input.nano} {input.pac}
        """

### Add an intermediate rule here with ragtag?
rule missing_gfa_to_fasta:
    input:
        prima = rules.hifiasmUL_reference_roi.output['prima'],
        hap1 = rules.hifiasmUL_reference_roi.output['hap1'],
        hap2 = rules.hifiasmUL_reference_roi.output['hap2'],
    output:
        prima = "map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_{scaffold}_primary.fa",
        hap1 = "map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_{scaffold}_hap1.fa",
        hap2 = "map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_{scaffold}_hap2.fa",
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input.prima} > {output.prima}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap1} > {output.hap1}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap2} > {output.hap2}
        """

rule missing_reassembled_scaffold:
    input:
        fa = "map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_{scaffold}_{assemb}.fa",
        ref = rules.missing_reference_chromosomes_extract.output['fa'],
    output:
        temp("map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_{scaffold}_{assemb}_scaffolded.fa")
    threads: 10
    params:
        "map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_{scaffold}_{assemb}/"
    conda:
        "../envs/ragtag.yaml"
    shell:
        """
        ragtag.py scaffold -t {threads} -o {params} -u -r -g 1 -m 2000000000 {input.ref} {input.fa}
        sed 's/RagTag/{wildcards.sample}_{wildcards.assemb}/g' {params}"ragtag.scaffold.fasta" > {output}
        rm -r {params}
        """

rule combine_scaffolded_missing_assemblies:
    input:
        expand("map_temp/{{species}}/{{sample}}/{{region}}/re-assemble/{{sample}}_{{hap}}_{scaffold}_{assemb}_scaffolded.fa", assemb = ['primary', 'hap1', 'hap2'], scaffold = SCAFFOLDS)
    output:
        "map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_roi_reassembled.fa",
    shell:
        """
        cat {input} > {output}
        """

use rule mapping_roi_flanking_genes_to_reference_genome as missing_reassembled_mapping_to_roi with:
    input:
        fa = rules.combine_scaffolded_missing_assemblies.output,
        flanks = rules.flanking_genes_in_region.output,
    output:
        sam = temp("map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_roi_reassembled.sam"),
        roi_region = "map_temp/{species}/{sample}/{region}/{sample}_{hap}_roi_reassembled_flanking-genes.txt"

rule missing_roi_sam_to_bed:
    input:
        sam = rules.missing_reassembled_mapping_to_roi.output['sam']
    output:
        bam = temp("map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_roi_reassembled.bam"),
        sort = temp("map_temp/{species}/{sample}/{region}/re-assemble/{sample}_{hap}_roi_reassembled_sort.bam"),
        bed = "results/{sample}/annotation/{region}/flanking/{sample}_{hap}_{species}_flanking_genes_reassembled_{region}.bed"
    threads: 2
    conda:
        "../envs/sambed.yaml"
    shell:
        """
        samtools view -@ {threads} -bh {input.sam} > {output.bam}
        samtools sort -@ {threads} -o {output.sort} {output.bam}
        bedtools bamtobed -i {output.sort} > {output.bed}
        """

## When missing or fragmented scaffolds found use reference to
## reassemble and redo the analysis 
rule check_missing_reassemble_flanking_genes:
    input:
        fa = rules.combine_scaffolded_missing_assemblies.output,
        bed = rules.missing_roi_sam_to_bed.output['bed'],
        info = "regions/{species}/{region}.info"
    output:
        "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}_missing.fa",
    conda:
        "../envs/sambed.yaml"
    params:
        outdir = "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}",
        fai = "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}_missing.fa.fai"
    shell:
        """
        python scripts/process_reassembled.py -i {input.info} -b {input.bed} -f {input.fa} -o {params.outdir}
        if [[ -f {params.fai} ]]; then
            rm {params.fai}
        fi
        """

rule check_fragmented_reassemble_flanking_genes:
    input:
        fa = rules.combine_scaffolded_missing_assemblies.output,
        bed = rules.missing_roi_sam_to_bed.output['bed'],
        info = "regions/{species}/{region}.info"
    output:
        "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}_fragmented.fa",
    conda:
        "../envs/sambed.yaml"
    params:
        "results/{sample}/annotation/{region}/{sample}_{hap}_{species}_{region}"
    shell:
        """
        python scripts/process_reassembled.py -i {input.info} -b {input.bed} -f {input.fa} -o {params}
        """

rule intact_analysis_done:
    input:
        "results/{sample}/annotation/{region}/single_libs/{sample}_{hap}_{species}_{region}_intacted_final_{i}.csv"
    output:
        "results/{sample}/annotation/{region}/intact/{sample}_{hap}_{species}_{region}_{i}.done"
    shell:
        """
        echo {wildcards.sample} {wildcards.hap} {wildcards.region} {wildcards.i} intacted > {output}
        """


## Annotation fragmented

rule fragmented_analysis_done:
    input:
        "results/{sample}/annotation/{region}/single_libs/{sample}_{hap}_{species}_{region}_fragmented_final_{i}.csv",
    output:
        "results/{sample}/annotation/{region}/fragmented/{sample}_{hap}_{species}_{region}_{i}.done"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools getfasta -fi {input.fa} -bed {input.bed} -s > {output}
        """

rule missing_analysis_done:
    input:
        "results/{sample}/annotation/{region}/single_libs/{sample}_{hap}_{species}_{region}_missing_final_{i}.csv"
    output:
        "results/{sample}/annotation/{region}/missing/{sample}_{hap}_{species}_{region}_{i}.done"
    shell:
        """
        echo {wildcards.sample} {wildcards.hap} {wildcards.region} {wildcards.i} missing > {output}
        """

def check_gaps(wildcards):
    with checkpoints.check_flanking_genes.get(**wildcards).output[0].open() as f:
        for line in f:
            col1, col2, col3, col4 = line.strip().split("\t")
            if col2 == "closed":
                return "results/{sample}/annotation/{region}/intact/{sample}_{hap}_{species}_{region}_{i}.done"
            elif col2 == "fragmented":
                return "results/{sample}/annotation/{region}/fragmented/{sample}_{hap}_{species}_{region}_{i}.done"
            else:
                return "results/{sample}/annotation/{region}/missing/{sample}_{hap}_{species}_{region}_{i}.done"

rule check_region:
    input:
        check_gaps
    output:
        temp("results/{sample}/annotation/{region}/flanking/{sample}_{hap}_{species}_{region}_{i}.annotated")
    shell:
        "cat {input} > {output}"


def combine_libraries(wildcards):
    checkpoint_output = checkpoints.libraries_format.get(**wildcards).output[0]
    return expand("results/{sample}/annotation/{region}/flanking/{sample}_{hap}_{species}_{region}_{i}.annotated",
        species=wildcards.species,
        region=wildcards.region,
        sample=wildcards.sample,
        hap=wildcards.hap,
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i)


rule multiple_libs_combine:
    input:
        combine_libraries
    output:
        "results/{sample}/annotation/{region}/final_results/{sample}_{hap}_{species}_{region}_combined.tsv",
    shell:
        """
        IFS=' ' read -r -a files <<< "{input}"
        num=${{#files[@]}}
        echo "Found $num input info file(s):" "${{files[@]}}"

        # 2) Initialize variables
        cdna="None"
        gdna="None"
        protein="None"

        # 3) Loop over each input‐line file, read its single line (or multiple lines),
        #    and set the appropriate variable based on $ref.
        for f in "${{files[@]}}"; do
            # (Optional check: ensure file actually exists)
            if [ ! -f "$f" ]; then
                echo "ERROR: Input file not found: $f" >&2
                exit 1
            fi
            while read -r samp hap lib ref roi; do
                echo "$samp"
                # Build the find pattern
                pattern="${{samp}}_${{hap}}*${{lib}}*final_${{ref}}.csv"
                file_path=$(find ./results/{wildcards.sample}/annotation/ -type f -name "$pattern" -print -quit)

                case "$ref" in
                    cDNA)
                        cdna="$file_path"
                        ;;
                    gDNA)
                        gdna="$file_path"
                        ;;
                    protein)
                        protein="$file_path"
                        ;;
                    *)
                        echo "ERROR: Unexpected ref value: $ref" >&2
                        exit 1
                        ;;
                esac
            done < "$f"
        done
        echo "python scripts/cgp_combine.py -c $cdna -g $gdna -p $protein -o {output}"
        python scripts/cgp_combine.py -c $cdna -g $gdna -p $protein -o {output} || touch {output}
        """

rule multiple_libs_refining:
    input:
        rules.multiple_libs_combine.output
    output:
        "results/{sample}/annotation/{region}/final_results/{species}_{region}_{sample}_{hap}_refined.bed",
    shell:
        """
        python scripts/cgp_refine4.py -i {input} -o {output} || touch {output}
        """

## Rules to detect gaps?
rule roi_figures:
    input:
        rules.multiple_libs_refining.output
    conda:
        "../envs/dna_viewer.yaml"
    output:
        "results/{sample}/annotation/{region}/final_results/{species}_{region}_{sample}_{hap}_refined.svg"
    params:
        "results/{sample}/annotation/{region}/final_results/{species}_{region}_{sample}_{hap}_refined"
    shell:
        """
        python scripts/dna_viewer_table.py -i {input} -o {params}
        """

rule refine_completed:
    input:
        files = combine_libraries,
        refined = rules.multiple_libs_refining.output
    output:
        temp("map_temp/{species}/{sample}/{region}/{species}_{region}_{sample}_{hap}_refine.txt")
    shell:
        """
        cat {input.files} > {output}
        """

rule annotation_summary:
    input:
        lambda wildcards: expand(
            "map_temp/{species}/{sample}/{region}/{species}_{region}_{sample}_{hap}_refine.txt",
            species=SPECIES,
            region=REGIONS,
            sample=wildcards.sample,
            hap=['hap1','hap2'],
        )
    output:
        "results/{sample}/annotation/{sample}_summary.txt"
    shell:
        """
        cat {input} > {output}
        """


### Check what this does
rule region_intact_accurate_coords:
    input:
        filtr = "results/{sample}/annotation/{region}/{hap}/{sample}_{hap}_{species}_{ref}_{region}_filter.tsv",
        bed = "results/{sample}/annotation/{region}/{hap}/{sample}_{hap}_{species}_{ref}_{region}_flanking_genes.bed",
    output:
        "results/{sample}/annotation/{region}/{hap}/intact/{sample}_{hap}_{species}_{ref}_{region}_final.tsv"
    benchmark:
        "results/{sample}/benchmarks/07_{sample}_{hap}_{species}_{ref}_{region}_11_coords_correct.bench"
    shell:
        """
        python scripts/coord_correction.py -b {input.bed} -f {input.filtr} -o {output}
        """







rule haplo_scaffolds_merge:
    input:
        hap1 = "results/{sample}/scaffolds/{sample}_hybridPN_hap1_scaffolded.fa",
        hap2 = "results/{sample}/scaffolds/{sample}_hybridPN_hap2_scaffolded.fa",
    output:
        "results/{sample}/midway/liftoff/{sample}_hybridPN_2haps.fasta"
    shell:
        """
        cat {input} > {output}
        """

rule haplo_genome_liftOff:
    input:
        genome = ancient("references/annotation/{ref}.fna"),
        gff = ancient("references/annotation/{ref}.gff"),
        fa = rules.haplo_scaffolds_merge.output,
    output:
        gff = "results/{sample}/scaffolds/{sample}_hybridPN_2haps_{ref}.gff",
        unmap = "results/{sample}/scaffolds/{sample}_hybridPN_2haps_{ref}_liftoff_unmapped.txt",
        outdir = temp(directory("results/{sample}/midway/liftoff/{ref}_liftoff_hybridPN")),
    benchmark:
        "results/{sample}/benchmarks/05_{sample}_hybridPN_{ref}_liftoff.bench"
    params:
        "results/{sample}/midway/liftoff/{sample}_hybridPN_2haps.fasta.fai"
    conda:
        "../envs/liftoff.yaml"
    threads: 10
    shell:
        """
        liftoff -g {input.gff} -o {output.gff} -p {threads} -u {output.unmap} -dir {output.outdir} {input.fa} {input.genome}
        if [[ -f {params} ]]; then
            rm {params}
        fi
        """










'''




rule haplo_genome_gff_correct:
    input:
        "results/{sample}/scaffolds/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_{ref}.gff",
    output:
        gff = "results/{sample}/scaffolds/{sample}_hybridPN_denovo_ntlink_{ref}.gff",
    benchmark:
        "results/{sample}/benchmarks/05_{sample}_hybridPN_{ref}_liftoff_correct.bench"
    conda:
        "../envs/genometools.yaml"
    shell:
        """
        gt gff3 -sort -tidy -retainids {input} > {output}
        echo "Complete liftoff analysis of sample {wildcards.sample} scaffolds"
        """


'''





















rule final_quast:
    input:
        contig_quast = "references/{scaffold}_quast/report.tsv",
        scaffold_quast = "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_{scaffold}_quast/report.tsv"
    output:
        stats = "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_ragtag_{scaffold}_stats.csv"
    benchmark:
        "results/{sample}/benchmarks/06_{sample}_hybridPN_{scaffold}_stats_collect.bench"
    shell:
        """
        grep "Assembly" {input.contig_quast} | awk '{{print $2}}' > {output}
        grep "contigs (>= 0 bp)" {input.contig_quast} | awk -F"\t" '{{print $2}}' >> {output}
        grep "N50" {input.contig_quast} | awk '{{print $2}}' >> {output}
        grep "contigs (>= 0 bp)" {input.scaffold_quast} | awk -F"\t" '{{print $2}}' >> {output}
        grep "N50" {input.scaffold_quast} | awk '{{print $2}}' >> {output}
        """



























'''
        samtools idxstats {output.sorted} > {output.chr} 
        samtools depth  {output.sorted}  |  awk '{{sum+=$3}} END {{ print "Average = ",sum/NR}}' > {output.cov}
'''






































rule prepare_report:
    input:
        "scripts/final_report_blank2.Rmd"
    output:
        "results/{sample}/info/{sample}_report.Rmd"
    shell:
        """
        cp {input} {output}
        """

rule full_report:
    input:
        sc = "results/{sample}/info/{sample}_report.Rmd",
        stats = "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_stats.csv",
        raw_stats = "results/{sample}/info/{sample}_filtered_reads_stats.txt",
        annotate = "results/{sample}/annotation/{sample}_{species}_{ref}_final_annotation.finish",
        busco = "results/{sample}/info/compare_{sample}_{ref}_figure.png",
    output:
        html = "results/{sample}/info/{sample}_{ref}_{species}_report.html"
    conda:
        "../envs/report.yaml"
    params:
        "results/{sample}/info/"
    log:
        "results/{sample}/logs/{sample}_{ref}_{species}_report.log"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.sc}', output_dir = '{params}', output_file=paste0('{wildcards.sample}_{wildcards.ref}_{wildcards.species}','_report'), output_format='html_document',params=list(subject='{wildcards.sample}',reference='{wildcards.ref}',species='{wildcards.species}'))" 2> {log}
        echo "Final report generated"
        """






onsuccess:
    shell("""
        echo "Complete analysis"
    """)

































































checkpoint hifiasm_haps_settings_busco:
    output:
        directory("results/{sample}/midway_others/{sample}_{hifiasm_mode}_busco/"),
    conda:
        "../envs/busco.yaml"
    params:
        busco = config["busco"],       
    threads: 20
    shell:
        """
        busco -m genome -i {input} -o {output} -l {params} -c {threads} -f
        """

rule haps_merge:
    input:
        rules.gfa_to_fasta_haps.output['hap1'],
        rules.gfa_to_fasta_haps.output['hap2'],
    output:
        "results/{sample}/midway/{sample}_hybridPN_2haps.fa"
    shell:
        """
        cat {input} > {output}
        """

checkpoint hifiasm_haps_busco:
    input:
        rules.haps_merge.output
    output:
        directory("results/{sample}/midway/{sample}_hybridPN_busco/"),
    conda:
        "../envs/busco.yaml"
    params:
        busco = config["busco"],       
    threads: 20
    shell:
        """
        busco -m genome -i {input} -o {output} -l {params} -c {threads} -f
        """











































rule filtered_raw_stats:
    input:
        nanopore = "results/{sample}/info/{sample}_nanopore_stat_filtered.txt",
        pacbio = "results/{sample}/info/{sample}_pacbio_stat_filtered.txt",
    output:
        "results/{sample}/info/{sample}_filtered_reads_stats.txt"
    shell:
        """
        echo "file num_seqs sum_len min_len avg_len max_len" > {output}
        cat {input} | sed '/file/d' | awk '{{print $1,$4,$5,$6,$7,$8}}' | sed 's/.*\\///g;s/,//g' >> {output}
        """

rule ntlink_fasta:
    input:
        quast = "results/{sample}/info/{sample}_hybridPN_denovo_quast/report.tsv",
        fa = "results/{sample}/midway/{sample}_hybridPN_denovo.fa",
    output:
        temp("results/{sample}/midway/ntlink/{sample}_hybridPN_denovo.fa")
    shell:
        """
        cp {input.fa} {output}
        """

rule ntLink_grid_run:
    input:
        fa = rules.ntlink_fasta.output,
        nano = rules.combine_filter_nanopore_fastq.output,
        pac = rules.filter_pacbio_fastq.output,
    output:
        fa = "results/{sample}/midway/ntlink/{sample}_hybridPN_denovo.fa.k{kval}.w{wval}.z1000.ntLink.scaffolds.gap_fill.fa",
    conda:
        "../envs/ntlink.yaml"
    threads: 10
    benchmark:
        "results/{sample}/benchmarks/03_{sample}_hybridPN_denovo_ntlink_k{kval}_w{wval}.bench"
    log:
        "results/{sample}/logs/02_{sample}_hybridPN_denovo_ntlink_k{kval}.w{wval}.log"
    shell:
        """
        ntLink scaffold gap_fill target={input.fa} reads="{input.nano} {input.pac}" t={threads} sensitive=True overlap=True extra_clean k={wildcards.kval} w={wildcards.wval} a=2 2> {log}

        filled=$(find "results/{wildcards.sample}/midway/ntlink/" -name "*gap_fill.fa" | wc -l)
        if [[ $filled -eq 9 ]]; then
            find "results/{wildcards.sample}/midway/ntlink/" \\( -name "{wildcards.sample}*tsv" -o -name "{wildcards.sample}*trim*" -o -name "{wildcards.sample}*aby*" \\) -exec rm {{}} + 
        fi
        """

rule ntlink_least_contigs:
    input:
        expand("results/{{sample}}/midway/ntlink/{{sample}}_hybridPN_denovo.fa.k{kval}.w{wval}.z1000.ntLink.scaffolds.gap_fill.fa", kval = KVALS, wval = WVALS),
    output:
        "results/{sample}/midway/{sample}_hybridPN_ntlink_least_contigs.txt"
    shell:
        """
        grep -c ">" {input} | sort -t ':' -k2 > {output}
        """

checkpoint extract_least_contigs:
    input:
        rules.ntlink_least_contigs.output
    output:
        directory("results/{sample}/midway/{sample}_hybridPN_ntLink2")
    shell:
        """
        least_contigs=$(head -1 {input} | cut -d":" -f1)

        mkdir -p {output}
        cp $least_contigs {output}
        """

def ntlink_output(wildcards):
    checkpoint_output = checkpoints.extract_least_contigs.get(**wildcards).output[0]
    return expand("results/{sample}/midway/{sample}_hybridPN_ntLink2/{i}.fa",
               sample=wildcards.sample,
               i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fa")).i)

rule ntlink_primary_final:
    input:
        ntlink_output
    output:
        "results/{sample}/scaffolds/ntlink/{sample}_hybridPN_ntlink.fasta"
    shell:
        """
        cat {input} > {output}
        echo "The rule ntlink_final complete for sample {wildcards.sample}"
        """
