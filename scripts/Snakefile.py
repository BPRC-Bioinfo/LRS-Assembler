import os
import glob

configfile: "configs/macaque.yaml"

SPECIES=config['species'].replace(" ", "_")
REGIONS = config['region']

NSAMPLES = list(config["nanopore"].keys())
PSAMPLES = list(config["pacbio"].keys())

SAMPLES = list(set(NSAMPLES) & set(PSAMPLES))

print (SAMPLES)

REFERENCE = list(set(config["reference"].keys()))

print (REFERENCE)

wildcard_constraints:
    sample="|".join(SAMPLES),
    hap = "hap1|hap2",
    ref = "|".join(REFERENCE)

KVALS = [24,32,40]
WVALS = [100, 250, 500]


rule all:
    input:
        expand("results/{sample}/info/{sample}_{ref}_{species}_report.html", sample = SAMPLES, ref = REFERENCE, species = SPECIES),
#        expand("results/{sample}/scaffolds/{sample}_hybridPN_denovo_{hap}_{ref}.fa", species = SPECIES, region = REGIONS, sample = SAMPLES, hap = ["hap1","hap2"], ref = REFERENCE),
#        expand("results/{sample}/scaffolds/ragtag/{ref}/{hap}/{sample}_hybridPN_denovo_{hap}.fa", species = SPECIES, region = REGIONS, sample = SAMPLES, hap = ["hap1","hap2"], ref = REFERENCE)
#        expand("annotated_{species}_{region}_{sample}_{hap}_{ref}.txt", species = SPECIES, region = REGIONS, sample = SAMPLES, hap = ["hap1","hap2"], ref = REFERENCE),
#        expand("flanking_genes/{species}/{region}.fasta", species = SPECIES, region = REGIONS),
#        expand("flanking_{species}_{region}_{sample}_{hap}_{ref}.txt", species = SPECIES, region = REGIONS, sample = SAMPLES, hap = [1,2], ref = REFERENCE),


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
    retries: 3
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
        for i in {input}; do
            cat $i
        done | seqkit sana | seqkit seq -m 5000 | seqkit rmdup -s -o {output}
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
        bam2fastq -o fastq/{wildcards.sample}_pacbio_raw {input.pacbio_bam}
        echo "The rule combine_pacbio_fastq complete for sample {wildcards.sample}"
        """

rule filter_pacbio_fastq:
    input:
        "fastq/{sample}_pacbio_raw.fastq.gz"
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
        filtered = "results/{sample}/raws/{sample}_{seqMachine}_5000.fastq.gz"
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
    benchmark:
        "results/{sample}/benchmarks/02_{sample}_hybridPN_hifiasm_denovo_assembly.bench"
    shell:
        """
        hifiasm -o results/{wildcards.sample}/hifiasm/{wildcards.sample}_hybridPN -t {threads} --ul {input.nano} {input.pac} 2> {log}
        echo "The rule hifiasmUL_de_novo complete for sample {wildcards.sample}"
        """

rule gfa_to_fasta:
    input:
        prima = ancient("results/{sample}/hifiasm/{sample}_hybridPN.bp.p_ctg.gfa"),
        hap1 = ancient("results/{sample}/hifiasm/{sample}_hybridPN.bp.hap1.p_ctg.gfa"),
        hap2 = ancient("results/{sample}/hifiasm/{sample}_hybridPN.bp.hap2.p_ctg.gfa"),
    output:
        prima = "results/{sample}/intermediates/{sample}_hybridPN_denovo.fa",
        hap1 = "results/{sample}/intermediates/{sample}_hybridPN_denovo_hap1.fa",
        hap2 = "results/{sample}/intermediates/{sample}_hybridPN_denovo_hap2.fa",
    benchmark:
        "results/{sample}/benchmarks/02_{sample}_hybridPN_hifiasm_denovo_gfa_to_fasta.bench"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input.prima} | sed "s/ptg/{wildcards.sample}_hybridPN_hifiasmUL_/g" > {output.prima}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap1} | sed "s/h1tg/{wildcards.sample}_hybridPN_hifiasmUL_hap1_/g" > {output.hap1}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap2} | sed "s/h2tg/{wildcards.sample}_hybridPN_hifiasmUL_hap2_/g" > {output.hap2}
        """

rule contigs_quast:
    input:
        "results/{sample}/intermediates/{sample}_hybridPN_denovo.fa",
    output:
        "results/{sample}/info/{sample}_hybridPN_denovo_quast/report.tsv"
    conda:
        "../envs/quast.yaml"
    benchmark:
        "results/{sample}/benchmarks/03_{sample}_hybridPN_quast.bench"
    params:
        "results/{sample}/info/{sample}_hybridPN_denovo_quast/"
    shell:
        """
        quast.py {input} -o {params}
        echo "The rule contig_quast complete for sample {wildcards.sample}"
        """

rule ntlink_fasta:
    input:
        quast = "results/{sample}/info/{sample}_hybridPN_denovo_quast/report.tsv",
        fa = "results/{sample}/intermediates/{sample}_hybridPN_denovo.fa",
    output:
        temp("results/{sample}/intermediates/ntlink/{sample}_hybridPN_denovo.fa")
    shell:
        """
        cp {input.fa} {output}
        """

rule ntLink_grid_run:
    input:
        fa = "results/{sample}/intermediates/ntlink/{sample}_hybridPN_denovo.fa",
        nano = "results/{sample}/raws/{sample}_nanopore_5000.fastq.gz",
        pac = "results/{sample}/raws/{sample}_pacbio_5000.fastq.gz",
    output:
        fa = "results/{sample}/intermediates/ntlink/{sample}_hybridPN_denovo.fa.k{kval}.w{wval}.z1000.ntLink.scaffolds.gap_fill.fa",
    conda:
        "../envs/ntlink.yaml"
    threads: 10
    benchmark:
        "results/{sample}/benchmarks/04_{sample}_hybridPN_denovo_ntlink_k{kval}_w{wval}.bench"
    log:
        "results/{sample}/logs/02_{sample}_hybridPN_denovo_ntlink_k{kval}.w{wval}.log"
    shell:
        """
        ntLink scaffold gap_fill target={input.fa} reads="{input.nano} {input.pac}" t={threads} sensitive=True overlap=True extra_clean k={wildcards.kval} w={wildcards.wval} a=2 2> {log}

        filled=$(find "results/{wildcards.sample}/intermediates/ntlink/" -name "*gap_fill.fa" | wc -l)
        if [[ $filled -eq 9 ]]; then
            find "results/{wildcards.sample}/intermediates/ntlink/" \\( -name "{wildcards.sample}*tsv" -o -name "{wildcards.sample}*trim*" -o -name "{wildcards.sample}*aby*" \\) -exec rm {{}} + 
        fi
        """

rule ntlink_least_contigs:
    input:
        expand("results/{{sample}}/intermediates/ntlink/{{sample}}_hybridPN_denovo.fa.k{kval}.w{wval}.z1000.ntLink.scaffolds.gap_fill.fa", kval = KVALS, wval = WVALS),
    output:
        "results/{sample}/intermediates/{sample}_hybridPN_ntlink_least_contigs.txt"
    shell:
        """
        grep -c ">" {input} | sort -t ':' -k2 > {output}
        """

checkpoint extract_least_contigs:
    input:
        "results/{sample}/intermediates/{sample}_hybridPN_ntlink_least_contigs.txt"
    output:
        directory("results/{sample}/intermediates/{sample}_hybridPN_ntLink2")
    shell:
        """
        least_contigs=$(head -1 {input} | cut -d":" -f1)

        mkdir -p {output}
        cp $least_contigs {output}
        """

def ntlink_output(wildcards):
    checkpoint_output = checkpoints.extract_least_contigs.get(**wildcards).output[0]
    return expand("results/{sample}/intermediates/{sample}_hybridPN_ntLink2/{i}.fa",
               sample=wildcards.sample,
               i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fa")).i)

rule ntlink_final:
    input:
        ntlink_output
    output:
        "results/{sample}/scaffolds/ntlink/{sample}_hybridPN_denovo_ntlink.fasta"
    shell:
        """
        cat {input} > {output}
        echo "The rule ntlink_final complete for sample {wildcards.sample}"
        """

rule prepare_genome_reference:
    output:
        "reference/{ref}.fna",
    conda:
        "../envs/datasets.yaml"
    retries: 3
    params:
        genome = lambda wildcards: config["reference"][wildcards.ref]['genome'],
        ref = lambda wildcards: config["reference"][wildcards.ref]['accession_number'],
    shell:
        """
        if [ ! -f {params.genome} ]; then
            echo "Genome file not found locally, downloading from NCBI"
            datasets download genome accession {params.ref} --filename {wildcards.ref}.zip --include genome
            unzip {wildcards.ref}.zip -d "reference"
            mv reference/ncbi_dataset/data/{wildcards.ref}/{wildcards.ref}*.fna {output}
        else
            echo "Using local genome file {params.genome}"
            cp {params.genome} {output}
        fi
        if [ -f "reference/README.md" ]; then
            rm reference/README.md
        fi
        echo "Complete prepare_genome_reference"
        """

rule prepare_gff_reference:
    output:
        "reference/{ref}.gff",
    conda:
        "../envs/datasets.yaml"
    retries: 3
    params:
        ref = lambda wildcards: config["reference"][wildcards.ref]['accession_number'],
        gff = lambda wildcards: config["reference"][wildcards.ref]['gff']
    shell:
        """
        if [ ! -f {params.gff} ]; then
            echo "Gff file not found locally, downloading from NCBI"
            datasets download genome accession {params.ref} --filename {wildcards.ref}.zip --include gff3
            unzip {wildcards.ref}.zip -d "reference"
            mv reference/ncbi_dataset/data/{wildcards.ref}/genomic.gff {output}
        else
            echo "Using local genome file {params.gff}"
            cp {params.gff} {output}
        fi
        if [ -f "reference/README.md" ]; then
            rm reference/README.md
        fi
        echo "Complete prepare_gff_reference"
        """

rule prepare_chromosome_reference:
    output:
        "reference/{ref}.info",
    retries: 3
    conda:
        "../envs/datasets.yaml"
    params:
        ref = lambda wildcards: config["reference"][wildcards.ref]['accession_number'],
        info = lambda wildcards: config["reference"][wildcards.ref]['chr_info']
    shell:
        """
        if [ ! -f {params} ]; then
            echo "Info file not found locally, downloading from NCBI"
            datasets download genome accession {params.ref} --filename {wildcards.ref}.zip --include seq-report
            unzip {wildcards.ref}.zip -d "reference"
            python scripts/json_parse.py reference/ncbi_dataset/data/{wildcards.ref}/sequence_report.jsonl | sed '/accession/d;/MT/d' | awk '{{print $1",@_chr"$2}}' > {output} 
        else
            echo "Using local genome file {params.info}"
            cp {params.info} {output}
        fi
        if [ -f "reference/README.md" ]; then
            rm reference/README.md
        fi
        echo "Complete prepare_chromosome_reference"
        """

rule ragTag:
    input:
        ref = ancient("reference/{ref}.fna"),
        fa = "results/{sample}/scaffolds/ntlink/{sample}_hybridPN_denovo_ntlink.fasta"
    output:
        "results/{sample}/scaffolds/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}.fa"
    params:
        "results/{sample}/scaffolds/ragtag/{ref}/"
    threads: 10
    benchmark:
        "results/{sample}/benchmarks/05_{sample}_hybridPN_{ref}_ragtag.bench"
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
        echo "Scaffold complete for sample {wildcards.sample} using {wildcards.ref}"
        """

rule ragtag_haplo:
    input:
        ref = "reference/{ref}.fna",
        fa = ancient("results/{sample}/intermediates/{sample}_hybridPN_denovo_{hap}.fa"),
    output:
        "results/{sample}/scaffolds/ragtag/{ref}/{hap}/{sample}_hybridPN_denovo_{hap}.fa"
    params:
        "results/{sample}/scaffolds/ragtag/{ref}/{hap}/"
    threads: 10
    conda:
        "../envs/ragtag.yaml"
    benchmark:
        "results/{sample}/benchmarks/05_{sample}_hybridPN_{hap}_{ref}_ragtag.bench"
    shell:
        """
        if [[ {wildcards.ref} != "none" ]]; then
            ragtag.py scaffold -t {threads} -o {params} -u {input.ref} {input.fa}
            sed 's/_RagTag//g' {params}"ragtag.scaffold.fasta" > {output}
        else
            cp {input.fa} {output}
        fi
        """

rule accession_chromosome:
    input:
        "reference/{ref}.info"
    output:
        "results/{sample}/scaffolds/ragtag/{ref}/{sample}_{ref}_accs_chrs.txt"
    shell:
        """
        sed 's/@/{wildcards.sample}/g' {input} > {output}
        """

rule rename_scaffolds:
    input:
        acc = "results/{sample}/scaffolds/ragtag/{ref}/{sample}_{ref}_accs_chrs.txt",
        fa = "results/{sample}/scaffolds/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}.fa"
    output:
        fa = "results/{sample}/scaffolds/{sample}_hybridPN_denovo_ntlink_{ref}.fa",
        acc = temp("{sample}_hybridPN_{ref}.temp")
    benchmark:
        "results/{sample}/benchmarks/06_{sample}_hybridPN_{ref}_rename_scaffolds.bench"
    shell:
        """
        awk -F',' '{{print $1}}' {input.acc} > {output.acc} 
        bash scripts/rename_awk.sh {output.acc} {wildcards.sample} {input.fa} {output.fa}
        while IFS="," read -r acc chrs; do
            sed "s/${{acc}}/${{chrs}}/g" -i {output.fa}
        done < {input.acc}
        """

rule rename_haplotypes:
    input:
        acc = "results/{sample}/scaffolds/ragtag/{ref}/{sample}_{ref}_accs_chrs.txt",
        fa = "results/{sample}/scaffolds/ragtag/{ref}/{hap}/{sample}_hybridPN_denovo_{hap}.fa"
    output:
        fa = "results/{sample}/scaffolds/{sample}_hybridPN_denovo_{hap}_{ref}.fa",
        acc = "{sample}_hybridPN_{ref}_{hap}.temp"
    benchmark:
        "results/{sample}/benchmarks/06_{sample}_hibridPN_{ref}_{hap}_rename_scaffolds.bench"
    shell:
        """
        awk -F',' '{{print $1}}' {input.acc} > {output.acc}
        bash scripts/rename_awk.sh {output.acc} {wildcards.sample} {input.fa} {output.fa}
        while IFS="," read -r acc chrs; do
            sed "s/${{acc}}/${{chrs}}_{wildcards.hap}/g" -i {output.fa}
        done < {input.acc}
        """

rule scaffold_quast:
    input:
        "results/{sample}/scaffolds/{sample}_hybridPN_denovo_ntlink_{ref}.fa",
    output:
        "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_{ref}_quast/report.tsv"
    conda:
        "../envs/quast.yaml"
    benchmark:
        "results/{sample}/benchmarks/06_{sample}_hybridPN_{ref}_quast_final.bench"
    params:
        "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_{ref}_quast/"
    shell:
        """
        quast.py {input} -o {params}
        echo "Complete quast analysis of sample {wildcards.sample} scaffolds"
        """

rule scaffold_busco:
    input:
        "results/{sample}/scaffolds/{sample}_hybridPN_denovo_ntlink_{ref}.fa"
    output:
        "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_{ref}_busco/run_primates_odb10/short_summary.txt"
    conda:
        "../envs/busco.yaml"
    log:
        "results/{sample}/logs/03_{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_busco.log"
    benchmark:
        "results/{sample}/benchmarks/06_{sample}_hybridPN_{ref}_busco.bench"
    params:
        busco = config["busco"],
        dir = "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_{ref}_busco/"
    threads: 20
    shell:
        """
        busco -m genome -i {input} -o {params.dir} -l {params.busco} -c {threads} -f 2> {log}
        echo "Complete BUSCO analysis of sample {wildcards.sample} scaffolds"
        """

rule liftOff:
    input:
        genome = "reference/{ref}.fna",
        gff = "reference/{ref}.gff",
        fa = "results/{sample}/scaffolds/{sample}_hybridPN_denovo_ntlink_{ref}.fa",
        chrInfo = "results/{sample}/ragtag/{ref}/{sample}_{ref}_accs_chrs.txt"
    output:
        temp = temp("results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_{ref}.gff"),
        chrs = temp("results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_{ref}_chrs.txt"),
        unmap = "results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_{ref}_liftoff_unmapped.txt",
        dir = temp(directory("results/{sample}/ragtag/{ref}/liftoff_hybridPN")),
    benchmark:
        "results/{sample}/benchmarks/06_{sample}_hybridPN_{ref}_liftoff.bench"
    conda:
        "../envs/liftoff.yaml"
    threads: 10
    shell:
        """
        if [[ {wildcards.ref} != "none" ]]; then
            grep ">" {input.fa} | grep "chr" | sed 's/>//g' | grep -f - {input.chrInfo} > {output.chrs}
            liftoff -g {input.gff} -o {output.temp} -p {threads} -u {output.unmap} -dir {output.dir} -chroms {output.chrs} {input.fa} {input.genome}
        else
            touch {output.temp} {output.unmap}
            mkdir -p {output.dir}
        fi
        """

rule gff_correct:
    input:
        "results/{sample}/ragtag/{ref}/{sample}_hybridPN_denovo_ntlink_{ref}.gff",
    output:
        gff = "results/{sample}/scaffolds/{sample}_hybridPN_denovo_ntlink_{ref}.gff",
    benchmark:
        "results/{sample}/benchmarks/06_{sample}_hybridPN_{ref}_liftoff_correct.bench"
    conda:
        "../envs/genometools.yaml"
    shell:
        """
        gt gff3 -sort -tidy -retainids {input} > {output}
        echo "Complete liftoff analysis of sample {wildcards.sample} scaffolds"
        """

rule final_quast_busco:
    input:
        gff = "results/{sample}/scaffolds/{sample}_hybridPN_denovo_ntlink_{ref}.gff",
        contig_quast = "results/{sample}/info/{sample}_hybridPN_denovo_quast/report.tsv",
        busco = "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_{ref}_busco/run_primates_odb10/short_summary.txt",
        scaffold_quast = "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_{ref}_quast/report.tsv"
    output:
        stats = "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_stats.csv"
    benchmark:
        "results/{sample}/benchmarks/08_{sample}_hybridPN_{ref}_stats_collect.bench"
    shell:
        """
        busco_db=$(grep "Total BUSCO" {input.busco} | awk '{{print $1}}')
        grep "Assembly" {input.contig_quast} | awk '{{print $2}}' > {output}
        grep "Complete BUSCO"  {input.busco} | awk -v var="$busco_db" '{{print $1/var*100}}' >> {output} 
        grep "contigs (>= 0 bp)" {input.contig_quast} | awk -F"\t" '{{print $2}}' >> {output}
        grep "N50" {input.contig_quast} | awk '{{print $2}}' >> {output}
        grep "contigs (>= 0 bp)" {input.scaffold_quast} | awk -F"\t" '{{print $2}}' >> {output}
        grep "N50" {input.scaffold_quast} | awk '{{print $2}}' >> {output}
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

rule full_report:
    input:
        sc = "scripts/final_report_blank.Rmd",
        stats = "results/{sample}/info/{sample}_hybridPN_denovo_ntlink_ragtag_{ref}_stats.csv",
        raw_stats = "results/{sample}/info/{sample}_filtered_reads_stats.txt",
        annotate = "results/{sample}/info/annotated_{species}_{sample}_{ref}.txt",
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
        Rscript -e "rmarkdown::render('{input.sc}', output_dir = '{params}', output_file=paste0('{wildcards.sample}_{wildcards.ref}_{wildcards.species}','_report'), output_format='html_document',params=list(sample='{wildcards.sample}',reference='{wildcards.ref}'))" 2> {log}
        echo "Final report generated"
        """

include:"../annotation/annotation.py"

rule running_annotation:
    input:
        expand("regions/{{species}}/{{sample}}/{region}/flanking/{{sample}}_{hap}_{{ref}}_{region}_annotation.finish", region = REGIONS, hap = ["hap1","hap2"])
    output:
        "results/{sample}/info/annotated_{species}_{sample}_{ref}.txt"
    shell:
        """
        # Processing the intermediate output
        cat {input} > {output}
        """

rule main_task_test:
    input:
        "checking_{species}_{region}_{sample}_{hap}_{ref}.txt"
    output:
        touch("checked_{species}_{region}_{sample}_{hap}_{ref}.csv")


onsuccess:
    shell("""
        if [ -d "reference/ncbi_dataset/" ]; then
            rm -rf reference/ncbi_dataset/
        fi
        if [ -f "reference/md5sum.txt" ]; then
            rm reference/md5sum.txt
        fi
    """)
