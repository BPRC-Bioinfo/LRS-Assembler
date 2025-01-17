import pandas as pd
from glob import glob

configfile: "configs/run-config.yaml"
REF = config['reference']

# Identify input files
NSAMPLES,NANOPORE, = glob_wildcards("raws/{sample}/nanopore/{files}.fastq.gz")
PSAMPLES,PACBIO, = glob_wildcards("raws/{sample}/pacbio/{file}.bam")

SAMPLES = set(NSAMPLES + PSAMPLES)

def get_filenames(pattern):
    return [os.path.splitext(os.path.basename(f))[0] for f in glob(pattern)]

# Create a dictionary to hold the presence of files without path and extension
samplesdict = {
    sample: {
        'pacbio': get_filenames(f"raws/{sample}/pacbio/*.bam"),
        'nanopore': get_filenames(f"raws/{sample}/nanopore/*.fastq.gz")
    }
    for sample in SAMPLES
}

# Determine the method for each sample
def final_outcome(sample):
    pacbio_present = bool(samplesdict[sample]['pacbio'])
    nanopore_present = bool(samplesdict[sample]['nanopore'])
    if pacbio_present and nanopore_present:
        return "hybridPN"
    elif pacbio_present:
        return "pacbio"
    elif nanopore_present:
        return "nanopore"

sample_to_outcome = {sample: final_outcome(sample) for sample in SAMPLES}

# Convert dictionary to data frame
data_df = pd.DataFrame.from_dict(sample_to_outcome, orient='index', columns=['Assembly']).reset_index()
data_df.columns = ['Sample', 'Assembly']

# Create DataFrame from references list
df_references = pd.DataFrame({'Reference': REF})

# Create cartesian product of samples and references
df_final = pd.merge(data_df.assign(key=1), df_references.assign(key=1), on='key').drop('key', axis=1)

# Condition if no ref
print(df_final)

# Function to extract raw files
def extractFastq(dir,fastqType):
    return samplesdict[dir][fastqType]

wildcard_constraints:
    sample="|".join(SAMPLES),
    hap = "hap1|hap2"

KVALS = [24,32,40]
WVALS = [100, 250, 500]

rule all:
    input:
        expand("results/{sample}/info/{sample}_{method}_denovo_ntlink_ragtag_{ref}_stats.csv", zip, sample=df_final['Sample'], method=df_final['Assembly'], ref=df_final['Reference']),

rule refPrepare:
    output:
        zip = temp("{ref}.zip"),
        fa = "genomes/{ref}/{ref}.fna",
        gff = "genomes/{ref}/{ref}.gff",
        json = "genomes/{ref}/{ref}.info"
    conda:
        "../envs/datasets.yaml"
    params:
        rm = "genomes/README.md",
        json = "genomes/ncbi_dataset/data/assembly_data_report.jsonl",
        json2 = "genomes/ncbi_dataset/data/dataset_catalog.json",
        dir = directory("genomes/ncbi_dataset/data/{ref}"),
    benchmark:
        "pipe_benchmarks/00_{ref}_download.bench"
    shell:
        """
        if [[ {wildcards.ref} != "none" ]]; then
            datasets download genome accession {wildcards.ref} --filename {output.zip} --include gff3,genome,seq-report
            unzip {output.zip} -d "genomes"
            # More than one fna?
            mv genomes/ncbi_dataset/data/{wildcards.ref}/{wildcards.ref}*.fna {output.fa}
            mv genomes/ncbi_dataset/data/{wildcards.ref}/genomic.gff {output.gff}
            mv genomes/ncbi_dataset/data/{wildcards.ref}/sequence_report.jsonl {output.json}
            rm -r {params}
        else
            mkdir -p {params.dir}
            touch {output}
            rm -r {params.dir}
        fi
        """

rule pacbio_index:
    input:
        bam = "raws/{sample}/pacbio/{pacbio}.bam",
    output:
        pbi = "raws/{sample}/pacbio/{pacbio}.bam.pbi",
    conda:
        "../envs/pbtk.yaml"
    benchmark:
        "pipe_benchmarks/01_{sample}_{pacbio}_index.bench"
    shell:
        """
        pbindex {input.bam}
        """

rule pacbio_bamToFastq:
    input:
        bam = "raws/{sample}/pacbio/{pacbio}.bam",
        pbi = "raws/{sample}/pacbio/{pacbio}.bam.pbi",
    output:
        temp("fastq/{sample}_{pacbio}.fastq.gz"),
    conda:
        "../envs/bam2fastx.yaml"
    benchmark:
        "pipe_benchmarks/01_{sample}_{pacbio}_bam2fastq.bench"
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
        for i in {input}; do
            echo $i
        done > {output}
        """

rule pacbio_combineFastq:
    input:
        files = lambda wildcards: expand("fastq/{sample}_{pacbio}.fastq.gz", sample = wildcards.sample, pacbio = extractFastq(wildcards.sample, "pacbio")),
        list = "results/{sample}/info/{sample}_pacbio_raws.txt"
    output:
        temp("results/{sample}/raws/{sample}_pacbio.fastq.gz")
    benchmark:
        "pipe_benchmarks/01_{sample}_pacbio_fastq_raw.bench"
    shell:
        """
        cat {input.files} > {output}
        """

rule nanopore_input:
    input:
        lambda wildcards: expand("raws/{sample}/nanopore/{nano}.gz", sample = wildcards.sample, nano = extractFastq(wildcards.sample, "nanopore")),
    output:
        "results/{sample}/info/{sample}_nanopore_raws.txt"
    benchmark:
        "pipe_benchmarks/01_{sample}_pacbio_fastq_input.bench"
    shell:
        """
        echo {input} | tr " " "\n" > {output}
        """

rule nanopore_combineFastq:
    input:
        files = lambda wildcards: expand("raws/{sample}/nanopore/{nano}.gz", sample = wildcards.sample, nano = extractFastq(wildcards.sample, "nanopore")),
        list = "results/{sample}/info/{sample}_nanopore_raws.txt"
    output:
        temp("results/{sample}/raws/{sample}_nanopore.fastq.gz")
    benchmark:
        "pipe_benchmarks/01_{sample}_nanopore_fastq_raw.bench"
    shell:
        """
        cat {input.files} > {output}
        """

rule filterReads:
    input:
        fq = "results/{sample}/raws/{sample}_{seqMachine}.fastq.gz",
    output:
        "results/{sample}/{sample}_{seqMachine}_5000.fastq.gz"
    wildcard_constraints:
        seqMachine="(nanopore|pacbio)"
    conda:
        "../envs/seqkit.yaml"
    benchmark:
        "pipe_benchmarks/01_{sample}_{seqMachine}_fastq_filter.bench"
    threads: 5
    shell:
        """
        # Filter with read length
        seqkit sana {input.fq} | seqkit seq -m 5000 | seqkit rmdup -s -o {output}
        """

rule readStatsRaw:
    input:
        raw = "results/{sample}/raws/{sample}_{seqMachine}.fastq.gz",
    output:
        raw = "results/{sample}/info/{sample}_{seqMachine}_stat_raw.txt",
    conda:
        "../envs/seqkit.yaml"
    benchmark:
        "pipe_benchmarks/01_{sample}_{seqMachine}_fastq_raw_stats.bench"
    threads: 5
    shell:
        """
        seqkit stats --threads {threads} {input.raw} > {output.raw}
        """

rule readStatsFiltered:
    input:
        filtered = "results/{sample}/{sample}_{seqMachine}_5000.fastq.gz"
    output:
        filtered = "results/{sample}/info/{sample}_{seqMachine}_stat_filtered.txt",
    wildcard_constraints:
        seqMachine="(nanopore|pacbio)"
    conda:
        "../envs/seqkit.yaml"
    benchmark:
        "pipe_benchmarks/01_{sample}_{seqMachine}_fastq_filter_stats.bench"
    threads: 5
    shell:
        """
        seqkit stats --threads {threads} {input.filtered} > {output.filtered}
        """

rule readStatsPacNano:
    input:
        pacraw = "results/{sample}/info/{sample}_nanopore_stat_raw.txt",
        pacfiltered = "results/{sample}/info/{sample}_nanopore_stat_filtered.txt",
        nanoraw = "results/{sample}/info/{sample}_pacbio_stat_raw.txt",
        nanofiltered = "results/{sample}/info/{sample}_pacbio_stat_filtered.txt",
    output:
        "results/{sample}/info/{sample}_hybridPN_5000_stat_filtered.txt"
    shell:
        """
        cat {input} > {output}
        """

rule hifiasmULDeNovo:
    input:
        nano = "results/{sample}/{sample}_nanopore_5000.fastq.gz",
        pac = "results/{sample}/{sample}_pacbio_5000.fastq.gz",
    output:
        prima = "results/{sample}/de_novo/hifiasm/{sample}_hybridPN.bp.p_ctg.gfa",
        hap1 = "results/{sample}/de_novo/hifiasm/{sample}_hybridPN.bp.hap1.p_ctg.gfa",
        hap2 = "results/{sample}/de_novo/hifiasm/{sample}_hybridPN.bp.hap2.p_ctg.gfa",
    conda:
        "../envs/hifiasm.yaml"
    log:
        "logs/01_{sample}_hybridPN_hifiasm_denovo.log"
    threads: 20
    benchmark:
        "pipe_benchmarks/02_{sample}_hybridPN_hifiasm_denovo_assembly.bench"
    shell:
        """
        hifiasm -o results/{wildcards.sample}/de_novo/hifiasm/{wildcards.sample}_hybridPN -t {threads} --ul {input.nano} {input.pac} 2> {log}
        """

rule hifiasmULDenovoGfaToFasta:
    input:
        prima = ancient("results/{sample}/de_novo/hifiasm/{sample}_hybridPN.bp.p_ctg.gfa"),
        hap1 = ancient("results/{sample}/de_novo/hifiasm/{sample}_hybridPN.bp.hap1.p_ctg.gfa"),
        hap2 = ancient("results/{sample}/de_novo/hifiasm/{sample}_hybridPN.bp.hap2.p_ctg.gfa"),
        filtered = "results/{sample}/info/{sample}_hybridPN_5000_stat_filtered.txt",
    output:
        prima = "results/{sample}/intermediates/{sample}_hybridPN_denovo.fa",
        hap1 = "results/{sample}/intermediates/{sample}_hybridPN_denovo_hap1.fa",
        hap2 = "results/{sample}/intermediates/{sample}_hybridPN_denovo_hap2.fa",
    benchmark:
        "pipe_benchmarks/02_{sample}_hybridPN_hifiasm_denovo_gfa_to_fasta.bench"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input.prima} | sed "s/ptg/{wildcards.sample}_hybridPN_hifiasmUL_denovo_/g" > {output.prima}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap1} | sed "s/h1tg/{wildcards.sample}_hybridPN_hifiasmUL_denovo_hap1_/g" > {output.hap1}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap2} | sed "s/h2tg/{wildcards.sample}_hybridPN_hifiasmUL_denovo_hap2_/g" > {output.hap2}
        """

rule hifiasmDeNovo:
    input:
        pac = "results/{sample}/{sample}_pacbio_5000.fastq.gz",
    output:
        prima = "results/{sample}/de_novo/hifiasm/{sample}_pacbio.bp.p_ctg.gfa",
        hap1 = "results/{sample}/de_novo/hifiasm/{sample}_pacbio.bp.hap1.p_ctg.gfa",
        hap2 = "results/{sample}/de_novo/hifiasm/{sample}_pacbio.bp.hap2.p_ctg.gfa",
    conda:
        "../envs/hifiasm.yaml"
    log:
        "logs/01_{sample}_pacbio_hifiasm_denovo.log"
    threads: 20
    benchmark:
        "pipe_benchmarks/02_{sample}_pacbio_hifiasm_denovo_assembly.bench"
    shell:
        """
        hifiasm -o results/{wildcards.sample}/de_novo/hifiasm/{wildcards.sample}_hybridPN -t {threads} {input.pac} 2> {log}
        """

rule hifiasmDenovoGfaToFasta:
    input:
        prima = ancient("results/{sample}/de_novo/hifiasm/{sample}_pacbio.bp.p_ctg.gfa"),
        hap1 = ancient("results/{sample}/de_novo/hifiasm/{sample}_pacbio.bp.hap1.p_ctg.gfa"),
        hap2 = ancient("results/{sample}/de_novo/hifiasm/{sample}_pacbio.bp.hap2.p_ctg.gfa"),
        filtered = "results/{sample}/info/{sample}_pacbio_stat_filtered.txt",
    output:
        prima = "results/{sample}/intermediates/{sample}_pacbio_denovo.fa",
        hap1 = "results/{sample}/intermediates/{sample}_pacbio_denovo_hap1.fa",
        hap2 = "results/{sample}/intermediates/{sample}_pacbio_denovo_hap2.fa",
    benchmark:
        "pipe_benchmarks/02_{sample}_pacbio_hifiasm_denovo_gfa_to_fasta.bench"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input.prima} | sed "s/ptg/{wildcards.sample}_pacbio_hifiasm_denovo_/g" > {output.prima}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap1} | sed "s/h1tg/{wildcards.sample}_pacbio_hifiasm_denovo_hap1_/g" > {output.hap1}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap2} | sed "s/h2tg/{wildcards.sample}_pacbio_hifiasm_denovo_hap2_/g" > {output.hap2}
        """

rule flyeDeNovo:
    input:
        "results/{sample}/{sample}_nanopore_5000.fastq.gz",
    output:
        "results/{sample}/de_novo/flye/{sample}_nanopore_flye/assembly.fasta"
    conda:
        "../envs/flye.yaml"
    log:
        "logs/01_{sample}_nanopore_denovo_flye.log"
    threads: 20
    params:
        "results/{sample}/de_novo/flye/{sample}_nanopore_flye"
    benchmark:
        "pipe_benchmarks/02_{sample}_nanopore_denovo_flye_assembly.bench"
    shell:
        """
        flye --nano-hq {input} --out-dir {params} -t {threads} 2> {log}
        """

rule flyeFasta:
    input:
        filtered = "results/{sample}/info/{sample}_nanopore_stat_filtered.txt",
        fa = "results/{sample}/de_novo/flye/{sample}_nanopore_flye/assembly.fasta"
    output:
        "results/{sample}/intermediates/{sample}_nanopore_denovo.fa",
    benchmark:
        "pipe_benchmarks/02_{sample}_nanopore_denovo_flye_fasta.bench"
    shell:
        """
        sed 's/contig/{wildcards.sample}_nanopore_flye_denovo/g' {input.fa} > {output}
        """

rule contigsQuast:
    input:
        "results/{sample}/intermediates/{sample}_{method}_denovo.fa",
    output:
        "results/{sample}/info/{sample}_{method}_denovo_quast/report.tsv"
    conda:
        "../envs/quast.yaml"
    benchmark:
        "pipe_benchmarks/03_{sample}_{method}_quast.bench"
    params:
        "results/{sample}/info/{sample}_{method}_denovo_quast/"
    shell:
        """
        quast.py {input} -o {params}
        """

rule ntlink_fasta:
    input:
        quast = "results/{sample}/info/{sample}_{method}_denovo_quast/report.tsv",
        fa = "results/{sample}/intermediates/{sample}_{method}_denovo.fa",
    output:
        temp("results/{sample}/intermediates/ntlink/{sample}_{method}_denovo.fa")
    shell:
        """
        cp {input.fa} {output}
        """

rule ntLinkPacnano:
    input:
        fa = "results/{sample}/intermediates/ntlink/{sample}_hybridPN_denovo.fa",
        nano = "results/{sample}/{sample}_nanopore_5000.fastq.gz",
        pac = "results/{sample}/{sample}_pacbio_5000.fastq.gz",
    output:
        fa = "results/{sample}/intermediates/ntlink/{sample}_hybridPN_denovo.fa.k{kval}.w{wval}.z1000.ntLink.scaffolds.gap_fill.fa",
    conda:
        "../envs/ntlink.yaml"
    threads: 10
    benchmark:
        "pipe_benchmarks/04_{sample}_hybridPN_denovo_ntlink_k{kval}_w{wval}.bench"
    log:
        "logs/02_{sample}_hybridPN_denovo_ntlink_k{kval}.w{wval}.log"
    shell:
        """
        ntLink scaffold gap_fill target={input.fa} reads="{input.nano} {input.pac}" t={threads} sensitive=True overlap=True extra_clean k={wildcards.kval} w={wildcards.wval} a=2 2> {log}

        filled=$(find "results/{wildcards.sample}/intermediates/ntlink/" -name "*gap_fill.fa" | wc -l)
        if [[ $filled -eq 9 ]]; then
            find "results/{wildcards.sample}/intermediates/ntlink/" \( -name "{wildcards.sample}*tsv" -o -name "{wildcards.sample}*trim*" -o -name "{wildcards.sample}*aby*" \) -exec rm {{}} + 
        fi
        """

rule ntLinkPabio:
    input:
        fa = "results/{sample}/intermediates/ntlink/{sample}_pacbio_denovo.fa",
        pac = "results/{sample}/{sample}_pacbio_5000.fastq.gz",
    output:
        fa = "results/{sample}/intermediates/ntlink/{sample}_pacbio_denovo.fa.k{kval}.w{wval}.z1000.ntLink.scaffolds.gap_fill.fa",
    conda:
        "../envs/ntlink.yaml"
    threads: 10
    benchmark:
        "pipe_benchmarks/04_{sample}_pacbio_denovo_ntlink_k{kval}_w{wval}.bench"
    log:
        "logs/02_{sample}_pacbio_denovo_ntlink_k{kval}.w{wval}.log"
    shell:
        """
        ntLink scaffold gap_fill target={input.fa} reads="{input.pac}" t={threads} sensitive=True overlap=True extra_clean k={wildcards.kval} w={wildcards.wval} a=2 2> {log}

        filled=$(find "results/{wildcards.sample}/intermediates/ntlink/" -name "*gap_fill.fa" | wc -l)
        if [[ $filled -eq 9 ]]; then
            find "results/{wildcards.sample}/intermediates/ntlink/" \( -name "{wildcards.sample}*tsv" -o -name "{wildcards.sample}*trim*" -o -name "{wildcards.sample}*aby*" \) -exec rm {{}} + 
        fi
        """

rule ntLinkNanopore:
    input:
        fa = "results/{sample}/intermediates/ntlink/{sample}_nanopore_denovo.fa",
        nano = "results/{sample}/{sample}_nanopore_5000.fastq.gz",
    output:
        fa = "results/{sample}/intermediates/ntlink/{sample}_nanopore_denovo.fa.k{kval}.w{wval}.z1000.ntLink.scaffolds.gap_fill.fa",
    conda:
        "../envs/ntlink.yaml"
    threads: 10
    benchmark:
        "pipe_benchmarks/04_{sample}_nanopore_denovo_ntlink_k{kval}_w{wval}.bench"
    log:
        "logs/02_{sample}_nanopore_denovo_ntlink_k{kval}.w{wval}.log"
    shell:
        """
        ntLink scaffold gap_fill target={input.fa} reads="{input.nano}" t={threads} sensitive=True overlap=True extra_clean k={wildcards.kval} w={wildcards.wval} a=2 2> {log}

        filled=$(find "results/{wildcards.sample}/intermediates/ntlink/" -name "*gap_fill.fa" | wc -l)
        if [[ $filled -eq 9 ]]; then
            find "results/{wildcards.sample}/intermediates/ntlink/" \( -name "{wildcards.sample}*tsv" -o -name "{wildcards.sample}*trim*" -o -name "{wildcards.sample}*aby*" \) -exec rm {{}} + 
        fi
        """

rule ntlink_least_contigs:
    input:
        expand("results/{{sample}}/intermediates/ntlink/{{sample}}_{{method}}_denovo.fa.k{kval}.w{wval}.z1000.ntLink.scaffolds.gap_fill.fa", kval = KVALS, wval = WVALS),
    output:
        "results/{sample}/intermediates/{sample}_{method}_ntlink_least_contigs.txt"
    shell:
        """
        grep -c ">" {input} | sort -t ':' -k2 > {output}
        """

checkpoint extract_least_contigs:
    input:
        "results/{sample}/intermediates/{sample}_{method}_ntlink_least_contigs.txt"
    output:
        directory("results/{sample}/intermediates/{sample}_{method}_ntLink2")
    shell:
        """
        least_contigs=$(head -1 {input} | cut -d":" -f1)

        mkdir -p {output}
        cp $least_contigs {output}
        """

def ntlink_output(wildcards):
    checkpoint_output = checkpoints.extract_least_contigs.get(**wildcards).output[0]
    return expand("results/{sample}/intermediates/{sample}_{method}_ntLink2/{i}.fa",
               sample=wildcards.sample,
               method=wildcards.method,
               i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fa")).i)

rule ntlink_final:
    input:
        ntlink_output
    output:
        "results/{sample}/linked/{sample}_{method}_denovo_ntlink.fasta"
    shell:
        """
        cat {input} > {output}
        """

rule ragTag:
    input:
        ref = "genomes/{ref}/{ref}.fna",
        fa = "results/{sample}/linked/{sample}_{method}_denovo_ntlink.fasta"
    output:
        "results/{sample}/ragtag/{ref}/{sample}_{method}_denovo_ntlink_ragtag_{ref}.fa"
    params:
        "results/{sample}/ragtag/{ref}/"
    threads: 10
    benchmark:
        "pipe_benchmarks/05_{sample}_{method}_{ref}_ragtag.bench"
    conda:
        "../envs/ragtag.yaml"
    shell:
        """
        if [[ {wildcards.ref} != "none" ]]; then
            ragtag.py scaffold -t {threads} -o {params} -u {input.ref} {input.fa}
            sed 's/_RagTag//g' {params}"ragtag.scaffold.fasta" > {output}
        else
            cp {input.fa} {output}
        fi
        """

rule ragTagHybrid:
    input:
        ref = "genomes/{ref}/{ref}.fna",
        fa = "results/{sample}/intermediates/{sample}_hybridPN_denovo_{hap}.fa",
    output:
        "results/{sample}/ragtag/{ref}/{hap}/{sample}_hybridPN_denovo_{hap}.fa"
    params:
        "results/{sample}/ragtag/{ref}/{hap}/"
    threads: 10
    conda:
        "../envs/ragtag.yaml"
    benchmark:
        "pipe_benchmarks/05_{sample}_hybridPN_{hap}_{ref}_ragtag.bench"
    shell:
        """
        if [[ {wildcards.ref} != "none" ]]; then
            ragtag.py scaffold -t {threads} -o {params} -u {input.ref} {input.fa}
            sed 's/_RagTag//g' {params}"ragtag.scaffold.fasta" > {output}
        else
            cp {input.fa} {output}
        fi
        """

rule acc_chr:
    input:
        "genomes/{ref}/{ref}.info"
    output:
        "results/{sample}/ragtag/{ref}/{sample}_{ref}_accs_chrs.txt"
    benchmark:
        "pipe_benchmarks/06_{sample}_{ref}_acc_chr.bench"
    shell:
        """
        python scripts/json_parse.py {input} | sed '/accession/d;/MT/d' | awk '{{print $1,"{wildcards.sample}_chr"$2}}' > {output}
        """

rule rename_scaffolds:
    input:
        acc = "results/{sample}/ragtag/{ref}/{sample}_{ref}_accs_chrs.txt",
        fa = "results/{sample}/ragtag/{ref}/{sample}_{method}_denovo_ntlink_ragtag_{ref}.fa"
    output:
        fa = "results/{sample}/contigs/{sample}_{method}_denovo_ntlink_{ref}.fa",
        acc = temp("{sample}_{method}_{ref}.temp")
    benchmark:
        "pipe_benchmarks/06_{sample}_{method}_{ref}_rename_scaffolds.bench"
    shell:
        """
        awk '{{print $1}}' {input.acc} > {output.acc} 
        bash scripts/rename_awk.sh {output.acc} {wildcards.sample} {input.fa} {output.fa}
        while read -r acc chrs; do
            sed "s/${{acc}}/${{chrs}}/g" -i {output.fa}
        done < {input.acc}
        """

rule rename_hapScaffolds:
    input:
        acc = "results/{sample}/ragtag/{ref}/{sample}_{ref}_accs_chrs.txt",
        fa = "results/{sample}/ragtag/{ref}/{hap}/{sample}_hybridPN_denovo_{hap}.fa"
    output:
        fa = "results/{sample}/contigs/{sample}_hybridPN_denovo_{hap}_{ref}.fa",
        acc = temp("{sample}_hybridPN_{ref}_{hap}.temp")
    benchmark:
        "pipe_benchmarks/06_{sample}_hibridPN_{ref}_{hap}_rename_scaffolds.bench"
    shell:
        """
        awk '{{print $1}}' {input.acc} | sed '1d' > {output.acc} 
        bash scripts/rename_awk.sh {output.acc} {wildcards.sample} {input.fa} {output.fa}
        while read -r acc chrs; do
            sed "s/${{acc}}/{wildcards.sample}_{wildcards.hap}_chr${{chrs}}/g" -i {output.fa}
            sed "s/contig/{wildcards.hap}_contig/g" -i {output.fa}
        done < {input.acc}
        """

rule contigsFinalQuast:
    input:
        "results/{sample}/contigs/{sample}_{method}_denovo_ntlink_{ref}.fa",
    output:
        "results/{sample}/info/{sample}_{method}_denovo_ntlink_{ref}_quast/report.tsv"
    conda:
        "../envs/quast.yaml"
    benchmark:
        "pipe_benchmarks/06_{sample}_{method}_{ref}_quast_final.bench"
    params:
        "results/{sample}/info/{sample}_{method}_denovo_ntlink_{ref}_quast/"
    shell:
        """
        quast.py {input} -o {params}
        """

rule buscoFinal:
    input:
        "results/{sample}/contigs/{sample}_{method}_denovo_ntlink_{ref}.fa"
    output:
        "results/{sample}/info/{sample}_{method}_denovo_ntlink_{ref}_busco/run_primates_odb10/short_summary.txt"
    conda:
        "../envs/busco.yaml"
    log:
        "logs/03_{sample}_{method}_denovo_ntlink_ragtag_{ref}_busco.log"
    benchmark:
        "pipe_benchmarks/06_{sample}_{method}_{ref}_busco.bench"
    params:
        "results/{sample}/info/{sample}_{method}_denovo_ntlink_{ref}_busco/"
    threads: 20
    shell:
        """
        busco -m genome -i {input} -o {params} -l primates_odb10 -c {threads} -f 2> {log}
        """

rule liftOff:
    input:
        genome = "genomes/{ref}/{ref}.fna",
        gff = "genomes/{ref}/{ref}.gff",
        fa = "results/{sample}/contigs/{sample}_{method}_denovo_ntlink_{ref}.fa",
        chrInfo = "results/{sample}/ragtag/{ref}/{sample}_{ref}_accs_chrs.txt"
    output:
        temp = temp("results/{sample}/ragtag/{ref}/{sample}_{method}_denovo_ntlink_{ref}.gff"),
        unmap = "results/{sample}/ragtag/{ref}/{sample}_{method}_denovo_ntlink_{ref}_liftoff_unmapped.txt",
        dir = temp(directory("results/{sample}/ragtag/{ref}/liftoff_{method}")),
    benchmark:
        "pipe_benchmarks/06_{sample}_{method}_{ref}_liftoff.bench"
    conda:
        "../envs/liftoff.yaml"
    threads: 15
    shell:
        """
        if [[ {wildcards.ref} != "none" ]]; then
            yChr=$(grep -c "chrY" {input.fa} || true)
            if [[ $yChr -eq 0 ]]; then
                sed '/chrY/d' {input.chrInfo} -i
            fi
            liftoff -g {input.gff} -o {output.temp} -p {threads} -u {output.unmap} -dir {output.dir} -chroms {input.chrInfo} {input.fa} {input.genome}
        else
            touch {output.temp} {output.unmap}
            mkdir -p {output.dir}
        fi
        """

rule gff_correct:
    input:
        "results/{sample}/ragtag/{ref}/{sample}_{method}_denovo_ntlink_{ref}.gff",
    output:
        gff = "results/{sample}/contigs/{sample}_{method}_denovo_ntlink_{ref}.gff",
    benchmark:
        "pipe_benchmarks/06_{sample}_{method}_{ref}_liftoff_correct.bench"
    conda:
        "../envs/genometools.yaml"
    shell:
        """
        gt gff3 -sort -tidy -retainids {input} > {output}
        """

rule assembly_coverage_pacbio:
    input:
        pac = "results/{sample}/{sample}_pacbio_5000.fastq.gz",
        fa = "results/{sample}/ragtag/{ref}/{sample}_pacbio_denovo_ntlink_ragtag_{ref}.fa"
    output:
        sam = temp("results/{sample}/ragtag/{ref}/{sample}_pacbio_denovo_ntlink_ragtag_{ref}_pac.sam"),
    conda:
        "../envs/minimap2.yaml"
    benchmark:
        "pipe_benchmarks/07_{sample}_pacbio_{ref}_assembly_coverage.bench"
    threads: 20
    shell:
        """
        minimap2 -ax map-hifi -t {threads} {input.fa} {input.pac} > {output}
        """

rule assembly_coverage_hybrid_pac:
    input:
        hap = expand("results/{{sample}}/contigs/{{sample}}_hybridPN_denovo_{hap}_{{ref}}.fa", hap = ["hap1","hap2"]),
        pac = "results/{sample}/{sample}_pacbio_5000.fastq.gz",
        fa = "results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}.fa"
    output:
        sam = temp("results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_pac.sam"),
    conda:
        "../envs/minimap2.yaml"
    benchmark:
        "pipe_benchmarks/07_{sample}_hybridPN_{ref}_assembly_coverage.bench"
    threads: 20
    shell:
        """
        minimap2 -ax map-hifi -t {threads} {input.fa} {input.pac} > {output}
        """

rule assembly_coverage_hybrid_nano:
    input:
        hap = expand("results/{{sample}}/contigs/{{sample}}_hybridPN_denovo_{hap}_{{ref}}.fa", hap = ["hap1","hap2"]),
        nano = "results/{sample}/{sample}_nanopore_5000.fastq.gz",
        fa = "results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}.fa"
    output:
        sam = temp("results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_nano.sam"),
    conda:
        "../envs/minimap2.yaml"
    benchmark:
        "pipe_benchmarks/07_{sample}_nanopore_{ref}_assembly_coverage.bench"
    threads: 20
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input.fa} {input.nano} > {output}
        """

rule hybrid_mapping_sam_bam:
    input:
        "results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_{pacORnano}.sam",
    output:
        temp("results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_{pacORnano}.bam"),
    conda:
        "../envs/minimap2.yaml"
    benchmark:
        "pipe_benchmarks/07_{sample}_hybridPN_{ref}_assembly_remap_{pacORnano}.bench"
    threads: 10
    shell:
        """
        samtools view -@ {threads} -bh {input} > {output}
        """

rule hybrid_mapping_sort:
    input:
        pac = "results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_{pacORnano}.bam",
    output:
        sort = temp("results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_{pacORnano}_sort.bam"),
        bai = temp("results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_{pacORnano}_sort.bam.bai"),
        cov = "results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_{pacORnano}_sort.coverage",
    conda:
        "../envs/minimap2.yaml"
    benchmark:
        "pipe_benchmarks/07_{sample}_hybridPN_{ref}_assembly_remap_{pacORnano}_sort.bench"
    threads: 10
    shell:
        """
        samtools sort -@ {threads} -o {output.sort} {input}
        samtools index {output.sort}
        samtools depth {output.sort} |  awk -v var={wildcards.pacORnano} '{{sum+=$3}} END {{ print "Average_"var" = "sum/NR }}' > {output.cov}
        """

rule hybrid_merge_pac_nano:
    input:
        bam = expand("results/{{sample}}/ragtag/{{ref}}/{{sample}}_hybridPN_denovo_ntlink_ragtag_{{ref}}_{pacORnano}_sort.bam", pacORnano = ["pac", "nano"]),
    output:
        temp("results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}.bam"),
    conda:
        "../envs/minimap2.yaml"
    benchmark:
        "pipe_benchmarks/07_{sample}_hybridPN_{ref}_assembly_remap_merge.bench"
    threads: 10
    shell:
        """
        samtools merge -@ {threads} -o {output} {input.bam}
        """

rule hybrid_mapping_coverage:
    input:
        expand("results/{{sample}}/ragtag/{{ref}}/{{sample}}_hybridPN_denovo_ntlink_ragtag_{{ref}}_{pacORnano}_sort.bam", pacORnano = ["pac", "nano"])
    output:
        "results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}.coverage",
    benchmark:
        "pipe_benchmarks/07_{sample}_hybridPN_{ref}_assembly_remap_coverage.bench"
    shell:
        """
        cat {input} | awk 'SUM+=$3;END{{print "Average_hybrid =", SUM}}' > {output}
        """

rule final_quast_busco:
    input:
        gff = "results/{sample}/contigs/{sample}_{method}_denovo_ntlink_{ref}.gff",
        quast = "results/{sample}/info/{sample}_{method}_denovo_quast/report.tsv",
        busco = "results/{sample}/info/{sample}_{method}_denovo_ntlink_{ref}_busco/run_primates_odb10/short_summary.txt",
        finalQuast = "results/{sample}/info/{sample}_{method}_denovo_ntlink_{ref}_quast/report.tsv"
    output:
        stats = "results/{sample}/info/{sample}_{method}_denovo_ntlink_ragtag_{ref}_stats.csv"
    benchmark:
        "pipe_benchmarks/08_{sample}_{method}_{ref}_stats_collect.bench"
    shell:
        """
        busco_db=$(grep "Total BUSCO" {input.busco} | awk '{{print $1}}')
        grep "Assembly" {input.quast} | awk '{{print $2}}' > {output}
        grep "Complete BUSCO"  {input.busco} | awk -v var="$busco_db" '{{print $1/var*100}}' >> {output} 
        grep "contigs (>= 0 bp)" {input.quast} | awk -F"\t" '{{print $2}}' >> {output}
        grep "N50" {input.quast} | awk '{{print $2}}' >> {output}
        grep "contigs (>= 0 bp)" {input.finalQuast} | awk -F"\t" '{{print $2}}' >> {output}
        grep "N50" {input.finalQuast} | awk '{{print $2}}' >> {output}
        """

