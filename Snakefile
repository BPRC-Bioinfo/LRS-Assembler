import pandas as pd

CHROMOSOMES = list(range(1,21)) + list("XY")
#CHROMOSOMES.append("Un")

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
    sample="|".join(SAMPLES),
    chr="|".join(map(str,CHROMOSOMES)),

rule all:
    input:
        expand("results/{sample}/de_novo/hifiasm/{sample}_pacnano_hifiasmUL.bp.p_ctg.gfa", sample = SAMPLES),
#        expand("results/{sample}/{sample}_{machine}_5000.fastq.gz", sample = SAMPLES, machine = MACHINES),
#        expand("results/{sample}/Assemblies/{sample}_{machine}_chrRef_{assembler}.fa", sample = SAMPLES, machine = MACHINES, assembler = ["hifiasm","flye"]),
#        expand("results/{sample}/Assemblies/{sample}_pacnano_chrRef_{assemblerUL}.fa", sample = SAMPLES, assemblerUL = ["hifiasmUL","flyeUL"]),
#        expand("results/{sample}/Assemblies/{sample}_{machine}_denovo_hifiasm.fa", sample = SAMPLES, machine = MACHINES),

#        expand("results/{sample}/de_novo/flye/{sample}_{machine}_flye/assembly.fasta", sample = SAMPLES, machine = MACHINES),

#        expand("results/{sample}/chr{chr}/{sample}_{machine}_chr{chr}.fastq.gz", sample = SAMPLES, machine = MACHINES, chr = CHROMOSOMES),
#        expand("results/{sample}/unmapped/{sample}_{machine}_unmapped.fastq.gz", sample = SAMPLES, machine = MACHINES),
#        expand("results/{sample}/assemblies_chrs/hifiasm/{sample}_{machine}_chr{chr}_hifiasm.fa", sample = SAMPLES, machine = MACHINES, chr = CHROMOSOMES),
        #expand("results/{sample}/assemblies_chrs/hifiasm/{sample}_pacnano_chr{chr}_hifiasm.fa", sample = SAMPLES, chr = CHROMOSOMES),
#        expand("results/{sample}/assemblies_chrs/flye/chr{chr}/{machine}/assembly.fasta", sample = SAMPLES, machine = MACHINES, chr = CHROMOSOMES),
#        expand("results/{sample}/assemblies_chrs/flye/chr{chr}/pacnano/assembly.fasta", sample = SAMPLES, chr = CHROMOSOMES),
#        expand("results/{sample}/Assemblies/{sample}_pacnano_chrRef_flye.fa", sample = SAMPLES),
#        expand("results/{sample}/assemblies_chrs/hifiasm/{sample}_{machine}_chr7_hifiasm.fa", sample = SAMPLES, machine = MACHINES, chr = CHROMOSOMES),
#        expand("results/{sample}/assemblies_chrs/hifiasm/{sample}_pacnano_chr7_hifiasm.fa", sample = SAMPLES, chr = CHROMOSOMES),
#        expand("results/{sample}/assemblies_chrs/flye/chr7/{machine}/assembly.fasta", sample = SAMPLES, machine = MACHINES, chr = CHROMOSOMES),
#        expand("results/{sample}/assemblies_chrs/flye/chr7/pacnano/assembly.fasta", sample = SAMPLES, chr = CHROMOSOMES),

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
        temp = temp("Mmul10/Mmul_chr{chr}.txt")
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        awk -F'\t' -v var={wildcards.chr} '$1=="chr"var {{print $5}}' {input.report} > {output.temp}
        awk -F'\t' -v var={wildcards.chr} '$1~"chr"var"_" {{print $5}}' {input.report} >> {output.temp}
        seqkit grep -f {output.temp} {input.fa} | sed "s/>/>chr{wildcards.chr}_/g" > {output.fa}
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
        #temp("results/{sample}/{sample}_{machine}_5000.fastq.gz")
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

rule unmappedReads:
    input:
        bam = "results/{sample}/{sample}_{machine}_sorted.bam",
        bai = "results/{sample}/{sample}_{machine}_sorted.bam.bai",
    output:
        bam = temp("results/{sample}/{sample}_{machine}_unmapped.bam"),
        bai = temp("results/{sample}/{sample}_{machine}_unmapped.bam.bai"),
        fq = "results/{sample}/unmapped/{sample}_{machine}_unmapped.fastq.gz",
    conda:
        "envs/minimap2.yaml"
    threads: 8
    shell:
        """
        samtools view -b -f 4 {input.bam} > {output.bam}
        samtools index {output.bam}
        samtools fastq -@ {threads} -0 {output.fq} {output.bam}
        """

rule bamInfo:
    input:
        sort = ancient("results/{sample}/{sample}_{machine}_sorted.bam"),
        bai = ancient("results/{sample}/{sample}_{machine}_sorted.bam.bai"),
    output:
        chr = "results/{sample}/info/{sample}_{machine}_readsMappedChrs.txt",
        cov = "results/{sample}/info/{sample}_{machine}_averageCoverage.txt",
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        samtools idxstats {input.sort} > {output.chr} 
        samtools depth {input.sort} | awk '{{sum+=$3}} END {{ print "{wildcards.sample} {wildcards.machine} genome_cov_avg = ",sum/NR}}' > {output.cov}
        """

checkpoint separateChrs:
    input:
        bam = "results/{sample}/{sample}_{machine}_sorted.bam",
        info = "results/{sample}/info/{sample}_{machine}_readsMappedChrs.txt",
        temp = "Mmul10/Mmul_chr{chr}.txt"
    output:
        directory("temp/{sample}_{machine}_{chr}")
    conda:
        "envs/minimap2.yaml"
    threads: 8
    benchmark:
        "benchmarks/05_{sample}_{machine}_chr{chr}.time"
    shell:
        """
        mkdir -p {output}
        for i in $(cat {input.temp}); do
            # Extract individual chrs
            samtools view -b {input.bam} "chr"{wildcards.chr}"_"$i > {output}/$i".bam"
            #nummap=$(samtools view -h {output}/$i".bam" | grep -v "@" | wc -l || true)
            #if [[ $nummap == "0" ]]; then
            #    rm {output}/$i".bam"
            #fi
        done
        """

def combineChrFragments(wildcards):
    outdir = checkpoints.separateChrs.get(**wildcards).output[0]
    bams = glob_wildcards(os.path.join(outdir, "{bam}.bam")).bam
    #print (bams)
    return expand(os.path.join(outdir, "{bam}.bam"),
        bam = bams)

rule rejoinChrs:
    input:
        combineChrFragments
    output:
        temp("{sample}_{machine}_chr{chr}.bam"),
    conda:
        "envs/minimap2.yaml"
    threads: 8
    shell:
        """
        ulimit -Sn 4096
        samtools merge -@ {threads} -o {output} {input} 
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
        "results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_{machine}_chr{chr}_hifiasm.bp.p_ctg.gfa",
        temp("results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_{machine}_chr{chr}_hifiasm.ec.bin"),
        temp("results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_{machine}_chr{chr}_hifiasm.ovlp.reverse.bin"),
        temp("results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_{machine}_chr{chr}_hifiasm.ovlp.source.bin"),
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

rule hifiasmChrGfaToFasta:
    input:
        ancient("results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_{machine}_chr{chr}_hifiasm.bp.p_ctg.gfa")
    output:
        "results/{sample}/assemblies_chrs/hifiasm/{sample}_{machine}_chr{chr}_hifiasm.fa"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} | sed "s/ptg/{wildcards.sample}_{wildcards.machine}_chr{wildcards.chr}_hifiasm_/g" > {output}
        """

rule hifiasmUlChr:
    input:
        nano = "results/{sample}/chr{chr}/{sample}_nanopore_chr{chr}.fastq.gz",
        pac = "results/{sample}/chr{chr}/{sample}_pacbio_chr{chr}.fastq.gz",
    output:
        "results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_pacnano_chr{chr}_hifiasmUL.bp.p_ctg.gfa"
    conda:
        "envs/hifiasm.yaml"
    log:
        "logs/07_{sample}_pacnano_chr{chr}_hifiasmUL.log"
    threads: 20
    benchmark:
        "benchmarks/07_{sample}_pacnano_chr{chr}_hifiasmUL.time"
    shell:
        """
        hifiasm -o results/{wildcards.sample}/assemblies_chrs/hifiasm/chr{wildcards.chr}/{wildcards.sample}_pacnano_chr{wildcards.chr}_hifiasmUL -t {threads} --ul {input.nano}  {input.pac} 2> {log}
        """

rule hifiasmULChrGfaToFasta:
    input:
        ancient("results/{sample}/assemblies_chrs/hifiasm/chr{chr}/{sample}_pacnano_chr{chr}_hifiasmUL.bp.p_ctg.gfa")
    output:
        "results/{sample}/assemblies_chrs/hifiasm/{sample}_pacnano_chr{chr}_hifiasm.fa"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} | sed "s/ptg/{wildcards.sample}_pacnano_chr{wildcards.chr}_hifiasmUL_/g" > {output}
        """

rule hifiasmDeNovo:
    input:
        "results/{sample}/{sample}_{machine}_5000.fastq.gz",
    output:
        "results/{sample}/de_novo/hifiasm/{sample}_{machine}_hifiasm.bp.p_ctg.gfa"
    conda:
        "envs/hifiasm.yaml"
    log:
        "logs/07_{sample}_{machine}_denovo_hifiasm.log"
    threads: 40
    benchmark:
        "benchmarks/07_{sample}_{machine}_denovo_hifiasm.time"
    shell:
        """
        hifiasm -o results/{wildcards.sample}/de_novo/hifiasm/{wildcards.sample}_{wildcards.machine}_hifiasm -t {threads} {input} 2> {log}
        """

rule hifiasmDenovoGfaToFasta:
    input:
        ancient("results/{sample}/de_novo/hifiasm/{sample}_{machine}_hifiasm.bp.p_ctg.gfa")
    output:
        "results/{sample}/Assemblies/{sample}_{machine}_denovo_hifiasm.fa"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} | sed "s/ptg/{wildcards.sample}_{wildcards.machine}_denovo_hifiasm_/g" > {output}
        """

rule hifiasmULDeNovo:
    input:
        nano = "results/{sample}/{sample}_nanopore_5000.fastq.gz",
        pac = "results/{sample}/{sample}_pacbio_5000.fastq.gz",
    output:
        "results/{sample}/de_novo/hifiasm/{sample}_pacnano_hifiasmUL.bp.p_ctg.gfa"
    conda:
        "envs/hifiasm.yaml"
    log:
        "logs/07_{sample}_pacnano_hifiasm_denovo.log"
    threads: 40
    benchmark:
        "benchmarks/07_{sample}_pacnano_hifiasm_denovo.time"
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
             flye --pacbio-hifi {input} --out-dir {params} -t {threads} 2> {log}
        elif [[ "{wildcards.machine}" == "nanopore" ]]; then
             flye --nano-hq {input} --out-dir {params} -t {threads} 2> {log}
        fi
        """

rule flyeUlChr:
    input:
        nano = "results/{sample}/chr{chr}/{sample}_nanopore_chr{chr}.fastq.gz",
        pac = "results/{sample}/chr{chr}/{sample}_pacbio_chr{chr}.fastq.gz",
    output:
        "results/{sample}/assemblies_chrs/flye/chr{chr}/pacnano/assembly.fasta"
    conda:
        "envs/flye.yaml"
    params:
        "results/{sample}/assemblies_chrs/flye/chr{chr}/pacnano"
    log:
        "logs/07_{sample}_pacnano_chr{chr}_flyeUL.log"
    threads: 40
    benchmark:
        "benchmarks/07_{sample}_pacnano_chr{chr}_flyeUL.time"
    shell:
        """
        #flye --nano-hq {input.pac} {input.nano} --iterations 0 --out-dir {output} --threads {threads} --scaffold --keep-haplotype
        flye --nano-hq {input.pac} {input.nano} --iterations 0 --out-dir {params} --threads {threads} 2> {log}
        flye --pacbio-hifi {input.pac} --resume-from polishing --out-dir {params} --threads {threads} 2>> {log}
        """

rule flyeChrRename:
    input:
        ancient("results/{sample}/assemblies_chrs/flye/chr{chr}/{machine}/assembly.fasta")
    output:
        "results/{sample}/assemblies_chrs/flye/{sample}_{machine}_chr{chr}_flye.fa"
    shell:
        """
        sed 's/contig/{wildcards.sample}_{wildcards.machine}_chr{wildcards.chr}_flye/g' {input} > {output}
        """

rule flyeUlChrRename:
    input:
        ancient("results/{sample}/assemblies_chrs/flye/chr{chr}/pacnano/assembly.fasta")
    output:
        "results/{sample}/assemblies_chrs/flye/{sample}_pacnano_chr{chr}_flyeUL.fa"
    shell:
        """
        sed 's/contig/{wildcards.sample}_pacnano_chr{wildcards.chr}_flyeUL/g' {input} > {output}
        """

rule flyeFullChrs:
    input:
        expand("results/{{sample}}/assemblies_chrs/flye/{{sample}}_{{machine}}_chr{chr}_flye.fa", chr = CHROMOSOMES)
    output:
        "results/{sample}/Assemblies/{sample}_{machine}_chrRef_flye.fa"
    shell:
        """
        cat {input} > {output}
        """

rule flyeULFullChrs:
    input:
        expand("results/{{sample}}/assemblies_chrs/flye/{{sample}}_pacnano_chr{chr}_flyeUL.fa", chr = CHROMOSOMES)
    output:
        "results/{sample}/Assemblies/{sample}_pacnano_chrRef_flyeUL.fa"
    shell:
        """
        cat {input} > {output}
        """

rule hifiasmFullChrs:
    input:
        expand("results/{{sample}}/assemblies_chrs/hifiasm/{{sample}}_{{machine}}_chr{chr}_hifiasm.fa", chr = CHROMOSOMES)
    output:
        "results/{sample}/Assemblies/{sample}_{machine}_chrRef_hifiasm.fa"
    shell:
        """
        cat {input} > {output}
        """

rule hifiasmULFullChrs:
    input:
        expand("results/{{sample}}/assemblies_chrs/hifiasm/{{sample}}_pacnano_chr{chr}_hifiasm.fa", chr = CHROMOSOMES)
    output:
        "results/{sample}/Assemblies/{sample}_pacnano_chrRef_hifiasmUL.fa"
    shell:
        """
        cat {input} > {output}
        """

rule flyeDeNovo:
    input:
        "results/{sample}/{sample}_{machine}_5000.fastq.gz",
    output:
        "results/{sample}/de_novo/flye/{sample}_{machine}_flye/assembly.fasta"
    conda:
        "envs/flye.yaml"
    log:
        "logs/07_{sample}_{machine}_denovo_flye.log"
    threads: 40
    params:
        "results/{sample}/de_novo/flye/{sample}_{machine}_flye"
    benchmark:
        "benchmarks/07_{sample}_{machine}_denovo_flye.time"
    shell:
        """
        if [[ "{wildcards.machine}" == "pacbio" ]]; then
             flye --pacbio-hifi {input} --out-dir {params} -t {threads} 2> {log}
        elif [[ "{wildcards.machine}" == "nanopore" ]]; then
             flye --nano-hq {input} --out-dir {params} -t {threads} 2> {log}
        fi
        """

rule flyeULDeNovo:
    input:
        nano = "results/{sample}/{sample}_nanopore_5000.fastq.gz",
        pac = "results/{sample}/{sample}_pacbio_5000.fastq.gz",
    output:
        "results/{sample}/de_novo/flye/{sample}_pacnano_flyeUL/assembly.fasta"
    conda:
        "envs/flye.yaml"
    params:
        "results/{sample}/de_novo/flye/{sample}_pacnano_flyeUL"
    log:
        "logs/07_{sample}_pacnano_flyeUL_denovo.log"
    threads: 40
    benchmark:
        "benchmarks/07_{sample}_pacnano_flyeUL_denovo.time"
    shell:
        """
        flye --nano-hq {input.pac} {input.nano} --iterations 0 --out-dir {params} --threads {threads} 2> {log}
        flye --pacbio-hifi {input.pac} --resume-from polishing --out-dir {params} --threads {threads} 2>> {log}
        """




rule evaluatePacbioChrAssembly:
    input:
        "results/{sample}/Assemblies/{sample}_pacnano_chrRef_flye.fa"
    output:
        busco = "results/{sample}/busco_{assembly}/short_summary.specific.primates_odb10.busco_{sample}.txt"
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
































rule evaluatePacbioChrAssembly2:
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

# do oncomplete remove temp files

## Flye settings
# scaffold = fill gaps with 100Ns
# --keep-haplotype breaks up contigs. Recover later stage
'''
