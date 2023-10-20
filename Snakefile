import pandas as pd

CHROMOSOMES = list(range(1,21)) + list("XY")
# Identify input files
NSAMPLES, NANOPORE, = glob_wildcards("raws/{sample}/nanopore/{files}.fastq.gz")
PSAMPLES, PACBIO ,= glob_wildcards("raws/{sample}/pacbio/{file}.bam")
SAMPLES = set(NSAMPLES + PSAMPLES)

samplesdict={}
# Loop through folder to identify normal and compressed files
for sample in SAMPLES:
    samplesdict[sample]={}
    BAM, = glob_wildcards("raws/"+ sample + "/pacbio/{bam}.bam")
    samplesdict[sample]['pacbio'] = BAM
    FASTQ, = glob_wildcards("raws/"+ sample + "/nanopore/{fastq}.fastq.gz")
    samplesdict[sample]['nanopore'] = FASTQ

# Function to extract fastq files
def extractFastq(dir,fastqType):
    return samplesdict[dir][fastqType]

# Identify sequencing machines
MACHINES = []
for key, inner_dict in samplesdict.items():
    for inner_key, inner_list in inner_dict.items():
        if inner_list:
            MACHINES.append(inner_key)

# Create new list
ASSEMBLIES = MACHINES.copy()
if "pacbio" in MACHINES and "nanopore" in MACHINES: 
    ASSEMBLIES.append("pacnano")

wildcard_constraints:
    sample="|".join(SAMPLES)

rule all:
    input:
        expand("results/{sample}/assemblies_chrs/hifiasm/{sample}_{machine}_chr{chr}_hifiasm.fa", sample = SAMPLES, machine = MACHINES, chr = CHROMOSOMES),
        expand("results/{sample}/assemblies_chrs/hifiasm/{sample}_pacnano_chr{chr}_hifiasm.fa", sample = SAMPLES, chr = CHROMOSOMES),
        expand("results/{sample}/de_novo/hifiasm/{sample}_pacnano_hifiasmUL.bp.p_utg.gfa", sample = SAMPLES)

rule refDownload:
    output:
        "Mmul10/GCA_003339765.3_Mmul_10_genomic.fna",
        "Mmul10/GCA_003339765.3_Mmul_10_assembly_report.txt"
    params:
        dir = "Mmul10/",
        gz = "Mmul10/GCA_003339765.3_Mmul_10_genomic.fna.gz"
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/339/765/GCA_003339765.3_Mmul_10/GCA_003339765.3_Mmul_10_genomic.fna.gz -P {params.dir}
        gunzip {params.gz}
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/339/765/GCA_003339765.3_Mmul_10/GCA_003339765.3_Mmul_10_assembly_report.txt -P {params.dir}
        """

rule extractChromosomes:
    input:
        fa = ancient("Mmul10/GCA_003339765.3_Mmul_10_genomic.fna"),
        report = ancient("Mmul10/GCA_003339765.3_Mmul_10_assembly_report.txt")
    output:
        fa = temp("Mmul10/Mmul10_chr{chr}.fasta"),
        temp = temp("Mmul10/Mmul_chr{chr}")
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        awk -F'\t' -v var={wildcards.chr} '$1=="chr"var {{print $5}}' {input.report} > {output.temp}
        seqkit grep -f {output.temp} {input.fa} > {output.fa}
        sed "s/>/>chr{wildcards.chr} /g" {output.fa} -i
        """

rule Mmul10:
    input:
        expand("Mmul10/Mmul10_chr{chr}.fasta", chr = CHROMOSOMES)
    output:
        "Mmul10/Mmul10_chrs.fasta"
    shell:
        """
        cat {input} > {output}
        """

rule pacbio_bamToFastq:
    input:
        bam = "raws/{sample}/pacbio/{pacbio}.bam",
        pbi = "raws/{sample}/pacbio/{pacbio}.bam.pbi",
    output:
        temp("fastq/{sample}_{pacbio}.fastq.gz"),
    conda:
        "envs/bam2fastx.yaml"
    benchmark:
        "benchmarks/01_{sample}_{pacbio}_bam2fastq.time"
    shell:
        """
        bam2fastq -o fastq/{wildcards.sample}_{wildcards.pacbio} {input.bam}
        """

rule pacbio_input:
    input:
        lambda wildcards: expand("fastq/{sample}_{pacbio}.fastq.gz", sample = wildcards.sample, pacbio = extractFastq(wildcards.sample, "pacbio")),
    output:
        "results/{sample}/info/{sample}_pacbio_raws.txt"
    shell:
        """
        echo {input} | tr " " "\n" > {output}
        """

rule pacbio_combineFastq:
    input:
        files = lambda wildcards: expand("fastq/{sample}_{pacbio}.fastq.gz", sample = wildcards.sample, pacbio = extractFastq(wildcards.sample, "pacbio")),
        list = "results/{sample}/info/{sample}_pacbio_raws.txt"
    output:
        "results/{sample}/raws/{sample}_pacbio.fastq.gz"
    shell:
        """
        cat {input.files} > {output}
        """

rule nanopore_input:
    input:
        lambda wildcards: expand("raws/{sample}/nanopore/{nano}.fastq.gz", sample = wildcards.sample, nano = extractFastq(wildcards.sample, "nanopore")),
    output:
        "results/{sample}/info/{sample}_nanopore_raws.txt"
    shell:
        """
        echo {input} | tr " " "\n" > {output}
        """

rule nanopore_fastq:
    input:
        files = lambda wildcards: expand("raws/{sample}/nanopore/{nano}.fastq.gz", sample = wildcards.sample, nano = extractFastq(wildcards.sample, "nanopore")),
        list = "results/{sample}/info/{sample}_nanopore_raws.txt"
    output:
        "results/{sample}/raws/{sample}_nanopore.fastq.gz"
    shell:
        """
        cat {input.files} > {output}
        """

rule rawReadStats:
    input:
        "results/{sample}/raws/{sample}_{machine}.fastq.gz",
    output:
        "results/{sample}/info/{sample}_{machine}_stat_raw.txt"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit stat {input} > {output}
        """

rule shortReadsFilter:
    input:
        stat = "results/{sample}/info/{sample}_{machine}_stat_raw.txt",
        fq = "results/{sample}/raws/{sample}_{machine}.fastq.gz",
    output:
        "results/{sample}/{sample}_{machine}_5000.fastq.gz"
    conda:
        "envs/seqkit.yaml"
    benchmark:
        "benchmarks/02_{sample}_{machine}_readsFilter.time"
    shell:
        """
        # Filter with read length
        seqkit seq -m 5000 {input.fq} | seqkit rmdup -s -o {output}
        """

rule filteredReadStats:
    input:
        "results/{sample}/{sample}_{machine}_5000.fastq.gz"
    output:
        "results/{sample}/info/{sample}_{machine}_5000_stat_filtered.txt"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit stat {input} > {output}
        """

rule refMapping:
    input:
        info = "results/{sample}/info/{sample}_{machine}_5000_stat_filtered.txt",
        fq = "results/{sample}/{sample}_{machine}_5000.fastq.gz",
        ref = "Mmul10/Mmul10_chrs.fasta"
    output:
        temp("{sample}_{machine}_Mmul10.sam")
    conda:
        "envs/minimap2.yaml"
    threads: 20
    log:
        "logs/03_{sample}_{machine}_mapMmul10.log"
    benchmark:
        "benchmarks/03_{sample}_{machine}_Mmul10.time"
    shell:
        """
        if [[ "{wildcards.machine}" == "pacbio" ]]; then
            minimap2 -ax map-hifi -t {threads} {input.ref} {input.fq} > {output} 2> {log}
        elif [[ "{wildcards.machine}" == "nanopore" ]]; then
            minimap2 -ax map-ont -t {threads} {input.ref} {input.fq} > {output} 2> {log}
        fi
        """

rule bamSort:
    input:
        "{sample}_{machine}_Mmul10.sam"
    output:
        bam = temp("{sample}_{machine}_Mmul10.bam"),
        sort = "results/{sample}/{sample}_{machine}_sorted.bam",
        bai = "results/{sample}/{sample}_{machine}_sorted.bam.bai",
    conda:
        "envs/minimap2.yaml"
    threads: 8
    benchmark:
        "benchmarks/04_{sample}_{machine}_mapMmul10_sort.time"
    shell:
        """
        # -b output bam -h include header -S sam_input -o output
        samtools view -@ {threads} -bh {input} > {output.bam}
        samtools sort -@ {threads} -o {output.sort} {output.bam}
        samtools index {output.sort}
        """

rule bamInfo:
    input:
        "results/{sample}/{sample}_{machine}_sorted.bam",
    output:
        chr = "results/{sample}/info/{sample}_{machine}_readsMappedChrs.txt",
        cov = "results/{sample}/info/{sample}_{machine}_averageCoverage.txt",
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        samtools idxstats {input} > {output.chr} 
        samtools depth {input} | awk '{{sum+=$3}} END {{ print "{wildcards.sample} {wildcards.machine} genome_cov_avg = ",sum/NR}}' > {output.cov}
        """

rule separateChrs:
    input:
        bam = "results/{sample}/{sample}_{machine}_sorted.bam",
        info = "results/{sample}/info/{sample}_{machine}_readsMappedChrs.txt",
    output:
        temp("{sample}_{machine}_chr{chr}.bam"),
    conda:
        "envs/minimap2.yaml"
    threads: 8
    benchmark:
        "benchmarks/05_{sample}_{machine}_chr{chr}.time"
    shell:
        """
        # Extract individual chrs
        samtools view -b {input.bam} chr{wildcards.chr} > {output}
        """

rule bamInfoChrs:
    input:
        "{sample}_{machine}_chr{chr}.bam",
    output:
        temp("results/{sample}/info/{sample}_{machine}_chr{chr}_averageCoverage.txt"),
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        samtools depth {input} | awk '{{sum+=$3}} END {{ print "{wildcards.sample} {wildcards.machine} chr{wildcards.chr}_cov_avg = ",sum/NR}}' > {output}
        """

rule averageCovChrs:
    input:
        expand("results/{{sample}}/info/{{sample}}_{{machine}}_chr{chr}_averageCoverage.txt", chr = CHROMOSOMES),
    output:
        "results/{sample}/info/{sample}_{machine}_chrsAverageCoverage.txt",
    shell:
        """
        cat {input} > {output}
        """

rule chrFastq:
    input:
        bam = "{sample}_{machine}_chr{chr}.bam",
        cov = "results/{sample}/info/{sample}_{machine}_chrsAverageCoverage.txt",
    output:
        "results/{sample}/chr{chr}/{sample}_{machine}_chr{chr}.fastq.gz",
    conda:
        "envs/minimap2.yaml"
    threads: 8
    benchmark:
        "benchmarks/06_{sample}_{machine}_chr{chr}_fastq.time"
    shell:
        """
        # Extract fastq
        samtools fastq -@ {threads} -0 {output} {input.bam}
        """

rule hifiasmChr:
    input:
        "results/{sample}/chr{chr}/{sample}_{machine}_chr{chr}.fastq.gz",
    output:
        "results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_{machine}_chr{chr}_hifiasm.bp.p_utg.gfa"
    conda:
        "envs/hifiasm.yaml"
    log:
        "logs/07_{sample}_{machine}_chr{chr}_hifiasm.log"
    threads: 20
    benchmark:
        "benchmarks/07_{sample}_{machine}_chr{chr}_hifiasm.time"
    shell:
        """
        hifiasm -o results/{wildcards.sample}/assemblies_chrs/hifiasm/chr{wildcards.chr}/{wildcards.sample}_{wildcards.machine}_chr{wildcards.chr}_hifiasm -t {threads} {input} 2> {log}
        """

rule hifiasmGfaToFasta:
    input:
        "results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_{machine}_chr{chr}_hifiasm.bp.p_utg.gfa"
    output:
        "results/{sample}/assemblies_chrs/hifiasm/{sample}_{machine}_chr{chr}_hifiasm.fa"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} | sed "s/utg/{wildcards.sample}_{wildcards.machine}_chr{wildcards.chr}_hifiasm/g" > {output}
        """

rule hifiasmUlChr:
    input:
        nano = "results/{sample}/chr{chr}/{sample}_nanopore_chr{chr}.fastq.gz",
        pac = "results/{sample}/chr{chr}/{sample}_pacbio_chr{chr}.fastq.gz",
    output:
        "results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_pacnano_chr{chr}_hifiasmUL.bp.p_utg.gfa"
    conda:
        "envs/hifiasm.yaml"
    log:
        "logs/07_{sample}_pacnano_chr{chr}_hifiasm.log"
    threads: 20
    benchmark:
        "benchmarks/07_{sample}_pacnano_chr{chr}_hifiasm.time"
    shell:
        """
        hifiasm -o results/{wildcards.sample}/assemblies_chrs/hifiasm/chr{wildcards.chr}/{wildcards.sample}_pacnano_chr{wildcards.chr}_hifiasmUL -t {threads} --ul {input.nano}  {input.pac} 2> {log}
        """

rule hifiasmULGfaToFasta:
    input:
        "results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_pacnano_chr{chr}_hifiasmUL.bp.p_utg.gfa"
    output:
        "results/{sample}/assemblies_chrs/hifiasm/{sample}_pacnano_chr{chr}_hifiasm.fa"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} | sed "s/utg/{wildcards.sample}_pacnano_chr{wildcards.chr}_hifiasm/g" > {output}
        """

rule hifiasmDeNovo:
    input:
        pac = "results/{sample}/{sample}_pacbio_5000.fastq.gz",
    output:
        "results/{sample}/de_novo/hifiasm/{sample}_hifiasm.bp.p_utg.gfa"
    conda:
        "envs/hifiasm.yaml"
    log:
        "logs/07_{sample}_hifiasm_denovo.log"
    threads: 40
    benchmark:
        "benchmarks/07_{sample}_hifiasm.time"
    shell:
        """
        hifiasm -o results/{wildcards.sample}/deNovo/hifiasm/{wildcards.sample}_hifiasm -t {threads} {input.pac} 2> {log}
        """

rule hifiasmULDeNovo:
    input:
        nano = "results/{sample}/{sample}_nanopore_5000.fastq.gz",
        pac = "results/{sample}/{sample}_pacbio_5000.fastq.gz",
    output:
        "results/{sample}/de_novo/hifiasm/{sample}_pacnano_hifiasmUL.bp.p_utg.gfa"
    conda:
        "envs/hifiasm.yaml"
    log:
        "logs/07_{sample}_pacnano_hifiasm_denovo.log"
    threads: 40
    benchmark:
        "benchmarks/07_{sample}_pacnano_hifiasm.time"
    shell:
        """
        hifiasm -o results/{wildcards.sample}/de_novo/hifiasm/{wildcards.sample}_pacnano_hifiasmUL -t {threads} --ul {input.nano} {input.pac} 2> {log}
        """















rule flyeChr:
    input:
        "results/{sample}/chr{chr}/{sample}_{machine}_chr{chr}.fastq.gz",
    output:
        "results/{sample}/assemblies_chrs/flye/chr{chr}/{machine}/assembly.fasta"
    conda:
        "envs/flye.yaml"
    params:
        "results/{sample}/assemblies_chrs/flye/chr{chr}/{machine}"
    log:
        "logs/07_{sample}_{machine}_chr{chr}_flye.log"
    threads: 40
    benchmark:
        "benchmarks/07_{sample}_{machine}_chr{chr}_flye.time"
    shell:
        """
        if [[ "{wildcards.machine}" == "pacbio" ]]; then
             flye --pacbio-hifi {input} --out-dir {params} --scaffold -t {threads} 2> {log}
        elif [[ "{wildcards.machine}" == "nanopore" ]]; then
             flye --nano-hq {input} --out-dir {params} --scaffold -t {threads} 2> {log}
             #flye --pacbio-hifi {input} --out-dir {output} --scaffold -t {threads} --keep-haplotype 2> {log}
             #flye --nano-hq {input} --out-dir {output} --scaffold -t {threads} --keep-haplotype 2> {log}
        fi
        """

rule flyeUlChr:
    input:
        nano = "results/{sample}/chr{chr}/{sample}_nanopore_chr{chr}.fastq.gz",
        pac = "results/{sample}/chr{chr}/{sample}_pacbio_chr{chr}.fastq.gz",
    output:
        "results/{sample}/assemblies_chrs/flyeUL/chr{chr}/pacnano/assembly.fasta"
    conda:
        "envs/flye.yaml"
    params:
        "results/{sample}/assemblies_chrs/flyeUL/chr{chr}/pacnano"
    log:
        "logs/07_{sample}_pacnano_chr{chr}_flyeUL.log"
    threads: 40
    benchmark:
        "benchmarks/07_{sample}_pacnano_chr{chr}_flyeUL.time"
    shell:
        """
        #flye --nano-hq {input.pac} {input.nano} --iterations 0 --out-dir {output} --threads {threads} --scaffold --keep-haplotype
        flye --nano-hq {input.pac} {input.nano} --iterations 0 --out-dir {params} --threads {threads} --scaffold 
        flye --pacbio-hifi {input.pac} --resume-from polishing --out-dir {params} --threads {threads} --scaffold 
        """

rule flyeChrRename:
    input:
        "results/{sample}/assemblies_chrs/flye/chr{chr}/{machine}/assembly.fasta"
    output:
        "results/{sample}/assemblies_chrs/flye/{sample}_{machine}_chr{chr}_flye.fa"
    shell:
        """
        sed 's/contig/{wildcards.sample}_{wildcards.machine}_{wildcards.chr}_/g' {input} > {output}
        """

rule flyeUlChrRename:
    input:
        "results/{sample}/assemblies_chrs/flyeUL/chr{chr}/pacnano/assembly.fasta"
    output:
        "results/{sample}/assemblies_chrs/flyeUL/{sample}_pacnano_chr{chr}_flyeUL.fa"
    shell:
        """
        sed 's/contig/{wildcards.sample}_pacnano_{wildcards.chr}_/g' {input} > {output}
        """

rule flyeFullChrs:
    input:
        expand("results/{{sample}}/assemblies_chrs/flye/{{sample}}_{{machine}}_chr{chr}_flye.fa", chr = CHROMOSOMES)
    output:
        "results/{sample}/assemblies_chrs/flye/all_chromosomes/{sample}_{machine}_allChromosomes.fa"
    shell:
        """
        cat {input} > {output}
        """

rule flyeULFullChrs:
    input:
        expand("results/{{sample}}/assemblies_chrs/flyeUL/{{sample}}_pacnano_chr{chr}_flyeUL.fa", chr = CHROMOSOMES)
    output:
        "results/{sample}/assemblies_chrs/flye/all_chromosomes/{sample}_pacnano_allChromosomes.fa"
    shell:
        """
        cat {input} > {output}
        """

rule hifiasmFullChrs:
    input:
        expand("results/{{sample}}/assemblies_chrs/hifiasm/{{sample}}_{{machine}}_chr{chr}_hifiasm.fa", chr = CHROMOSOMES)
    output:
        "results/{sample}/assemblies_chrs/hifiasm/all_chromosomes/{sample}_{machine}_allChromosomes.fa"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} | sed "s/utg/{wildcards.sample}_{wildcards.machine}_chr{wildcards.chr}_/g" > {output}
        """































rule canuChr:
    input:
        "results/{sample}/chr{chr}/{sample}_{machine}_chr{chr}.fastq.gz",
    output:
        directory("results/{sample}/assemblies_chrs/canu/chr{chr}/{machine}")
    conda:
        "envs/canu.yaml"
    log:
        "logs/07_{sample}_{machine}_chr{chr}_canu.log"
    threads: 20
    benchmark:
        "benchmarks/07_{sample}_{machine}_chr{chr}_canu.time"
    shell:
        """
        if [[ "{wildcards.machine}" == "pacbio" ]]; then
             canu -p {wildcards.sample}_{wildcards.machine}_chr{chr}_canu -d {output} -pacbio-hifi {input} 2> {log}
        elif [[ "{wildcards.machine}" == "nanopore" ]]; then
             canu -p {wildcards.sample}_{wildcards.machine}_chr{chr}_canu -d {output} -nanopore {input} 2> {log}
        fi
        """

rule canuUlChr:
    input:
        nano = "results/{sample}/chr{chr}/{sample}_nanopore_chr{chr}.fastq.gz",
        pac = "results/{sample}/chr{chr}/{sample}_pacbio_chr{chr}.fastq.gz",
    output:
        directory("results/{sample}/assemblies_chrs/canuUL/chr{chr}")
    conda:
        "envs/flye.yaml"
    log:
        "logs/07_{sample}_pacnano_chr{chr}_canuUL.log"
    threads: 20
    benchmark:
        "benchmarks/07_{sample}_pacnano_chr{chr}_canuUL.time"
    shell:
        """
        canu -p {wildcards.sample}_{wildcards.machine}_chr{chr}_canuUl -d {output} -nanopore {input.nano} -pacbio-hifi {input.pac} 2> {log}
        """





























rule pacbioAllChrs:
    input:
        expand("results/{{sample}}/chr{chr}/{{sample}}_pacbio_chr{chr}.bp.p_utg.fa", chr = CHROMOSOMES),
    output:
        "results/{sample}/{sample}_pacbio_contigs.fasta"
    shell:
        """
        cat {input} > {output}
        """

rule pacnanoChrAssembly:
    input:
        pacbio = "results/{sample}/chr{chr}/{sample}_pacbio_chr{chr}.fastq.gz",
        nanopore = "results/{sample}/chr{chr}/{sample}_nanopore_chr{chr}.fastq.gz",
    output:
        "results/{sample}/chr{chr}/{sample}_pacnano_chr{chr}.bp.p_utg.gfa"
    conda:
        "envs/hifiasm.yaml"
    log:
        "logs/07_{sample}_pacnano_chr{chr}_assembly.log"
    threads: 20
    benchmark:
        "benchmarks/07_{sample}_pacnano_chr{chr}_assembly.time"
    shell:
        """
        hifiasm -o results/{wildcards.sample}/chr{wildcards.chr}/{wildcards.sample}_pacnano_chr{wildcards.chr} -t {threads} --ul {input.nanopore} {input.pacbio} 2> {log}
        """

rule pacnanoGfaToFasta:
    input:
        "results/{sample}/chr{chr}/{sample}_pacnano_chr{chr}.bp.p_utg.gfa"
    output:
        "results/{sample}/chr{chr}/{sample}_pacnano_chr{chr}.bp.p_utg.fa"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} | sed "s/utg/{wildcards.sample}_pacnano_chr{wildcards.chr}_/g" > {output}
        """

rule pacnanoAllChrs:
    input:
        expand("results/{{sample}}/chr{chr}/{{sample}}_pacnano_chr{chr}.bp.p_utg.fa", chr = CHROMOSOMES),
    output:
        "results/{sample}/{sample}_pacnano_contigs.fasta"
    shell:
        """
        cat {input} > {output}
        """

rule evaluatePacbioChrAssembly:
    input:
        "results/{sample}/{sample}_{assembly}_contigs.fasta"
    output:
        busco = "results/{sample}/busco_{assembly}/short_summary.specific.primates_odb10.busco_{assembly}.txt"
    conda:
        "envs/busco.yaml"
    params:
        "results/{sample}/busco_{assembly}/"
    log:
        "logs/08_{sample}_busco_{assembly}.log"
    threads: 20
    benchmark:
        "benchmarks/08_{sample}_busco_{assembly}.time"
    shell:
        """
        busco -m genome -i {input} -o {params} -l primates_odb10 -c {threads} -f 2> {log}
        """

'''
def final_out():
    if 'pacbio' in MACHINES and 'nanopore' in MACHINES:
        return(expand("results/{sample}/busco_pacnano/short_summary.specific.primates_odb10.busco_pacnano.txt", sample = SAMPLES))
    elif 'pacbio' in MACHINES:
        return(expand("results/{sample}/busco_pacbio/short_summary.specific.primates_odb10.busco_pacbio.txt", sample = SAMPLES))
    elif 'nanopore' in MACHINES:
        return(expand("results/{sample}/busco_pacbio/short_summary.specific.primates_odb10.busco_nanopore.txt", sample = SAMPLES))

rule all:
    input:
        final_out()
'''
