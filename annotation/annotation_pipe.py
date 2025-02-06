
## Flanking genes detection
rule download_flanking_genes:
    params:
        left = lambda wildcards: REGIONS[wildcards.region]['left_flank'],
        right = lambda wildcards: REGIONS[wildcards.region]['right_flank'],
        ref = config['species']
    output:
        temp = temp("flanking_genes/{species}/temp_{region}.fa"),
        final = "flanking_genes/{species}/{region}.fasta",
        info = "flanking_genes/{species}/info_{region}.txt"
    retries: 3
    conda:
        "../envs/flanking.yaml"
    shell:
        """
        left_flank={params.left}
        touch {output.info}

        if [[ ! -z $left_flank ]]; then
            datasets download gene symbol {params.left} --taxon "{params.ref}" --include gene --filename {params.left}.zip || (echo "Warning: Download left flanking failed! Check if gene or sepecies is correct.")
            unzip {params.left}.zip -d {params.left}
            cat {params.left}/ncbi_dataset/data/gene.fna >> {output.temp}
            rm -r {params.left} {params.left}.zip
            echo "left:" {params.left} >> {output.info}
            sleep 3
        fi

        right_flank={params.right}
        # Right flanking gene
        if [[ ! -z $right_flank ]]; then
            datasets download gene symbol {params.right} --taxon "{params.ref}" --include gene --filename {params.right}.zip || (echo "Warning: Download right flanking failed! Check if gene or sepecies is correct.")
            unzip {params.right}.zip -d {params.right}
            cat {params.right}/ncbi_dataset/data/gene.fna >> {output.temp}
            rm -r {params.right}/ {params.right}.zip
            echo "right:" {params.right} >> {output.info}
            sleep 2
        fi
        
        sed -E "s/ \\[.*//g;s/ /_/g" {output.temp} > {output.final}
        echo "Complete download flanking gene"
        """

rule locate_flanking_genes:
    input:
        fa = "results/{sample}/scaffolds/{sample}_hybridPN_denovo_{hap}_{ref}.fa",
        flanks = "flanking_genes/{species}/{region}.fasta"
    output:
        sam = "map_temp/{species}/{region}/{sample}_{hap}_{ref}_{region}_flank_genes.sam",
        flank_region = "map_temp/{species}/{region}/{sample}_{hap}_{ref}_{region}_flank_genes.txt"
    threads: 6
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        #minimap2 -ax splice -t {threads} {input.fa} {input.flanks} > {output.sam}
        minimap2 -ax asm5 -t {threads} {input.fa} {input.flanks} > {output.sam}
        # Filter low quality mapping flanking genes?
        # cut -f1-5 {output.sam} | grep -v "@" | awk -F'\t' '$5 > 0' > {output.flank_region}
        cut -f1-5 {output.sam} | grep -v "@" | awk '$3!="*"' > {output.flank_region}
        echo "Mapped flanking gene to {wildcards.sample} {wildcards.hap}"
        """

checkpoint check_flanking_genes:
    input:
        loc = "map_temp/{species}/{region}/{sample}_{hap}_{ref}_{region}_flank_genes.txt",
        info = "flanking_genes/{species}/info_{region}.txt"
    output:
        status = "regions/{species}/{sample}/{region}/{sample}_{hap}_{ref}_{region}_status.csv",
        tmp = temp("regions/{species}/{sample}/{region}/{sample}_{hap}_{ref}_{region}_status.tmp"),
    shell:
        """
        # Extract left and right flanking genes
        left_flank=$(grep 'left' {input.info} | awk -F': ' '{{print $2}}' || true)
        right_flank=$(grep 'right' {input.info} | awk -F': ' '{{print $2}}' || true)

        # Initialize variables to track the contigs for left and right flanking gene
        left_contigs=""
        right_contigs=""

        # Find contigs for left and right genes
        while IFS=$'\t' read -r region_info _ contig position _; do
            gene_name=$(echo $region_info | awk -F'_' '{{print $NF}}')

            if [ "$gene_name" = "$left_flank" ]; then
                left_contigs="$left_contigs $contig"
            elif [ "$gene_name" = "$right_flank" ]; then
                right_contigs="$right_contigs $contig"
            fi
        done < {input.loc}

        # Remove duplicates from the contig lists
        left_contigs=$(echo "$left_contigs" | tr ' ' '\n' | sort -u | tr '\n' ' ')
        right_contigs=$(echo "$right_contigs" | tr ' ' '\n' | sort -u | tr '\n' ' ')

        # Check if both flanking genes are found
        if [ -n "$left_flank" ] && [ -n "$right_flank" ]; then
            common_contigs=$(comm -12 <(echo "$left_contigs" | tr ' ' '\n' | grep -v '^$' | sort) <(echo "$right_contigs" | tr ' ' '\n' | grep -v '^$' | sort))
            # Find unique contigs to the left
            if [ -n "$common_contigs" ]; then
                for contig in $common_contigs; do
                    echo -e "flanks\tclosed\t$left_flank/$right_flank\t$contig" >> {output.tmp}
                done
            else
                left_unique=$(comm -23 <(echo "$left_contigs" | tr ' ' '\n' | grep -v '^$' | sort) <(echo "$right_contigs" | tr ' ' '\n' | grep -v '^$' | sort)| tr '\n' ' ')
                right_unique=$(comm -13 <(echo "$left_contigs" | tr ' ' '\n' | grep -v '^$' | sort) <(echo "$right_contigs" | tr ' ' '\n' | grep -v '^$' | sort)| tr '\n' ' ')

                if [[ -z "$left_unique" && -z "$right_unique" ]]; then
                    echo -e "flanks\tmissing\t*\t*" >> {output.tmp}
                else
                    echo -e "flanks\tfragmented\t$left_flank/$right_flank\t$left_unique/$right_unique" >> {output.tmp}
                fi
            fi
        elif [ -n "$left_flank" ]; then
            if [ -n "$left_contigs" ]; then

                for contig in $left_contigs; do
                    echo -e "telomere\tclosed\t$left_flank/\t$contig" >> {output.tmp}
                done    
            fi
        elif [ -n "$right_flank" ]; then
            if [ -n "$right_contigs" ]; then
                for contig in $left_contigs; do
                    echo -e "telomere\tclosed\t/$right_flank\t$contig" >> {output.tmp}
                done
            fi
        fi
        awk '!s[$0]++' {output.tmp} > {output.status}
        echo "Checking flanking gene to {wildcards.sample} {wildcards.hap}"
        """

## Processing complete fragments
rule region_sam_to_sort_bam:
    input:
        sam = "map_temp/{species}/{region}/{sample}_{hap}_{ref}_{region}_flank_genes.sam",
    output:
        bam = temp("map_temp/{species}/{region}/{sample}_{hap}_{ref}_{region}_flank_genes.bam"),
        sort = temp("map_temp/{species}/{region}/{sample}_{hap}_{ref}_{region}_flank_genes_sort.bam"),
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -@ {threads} -bh {input.sam} > {output.bam}
        samtools sort -@ {threads} -o {output.sort} {output.bam}
        """

rule region_intact_bed:
    input:
        sam = "map_temp/{species}/{region}/{sample}_{hap}_{ref}_{region}_flank_genes.sam",
        bam = "map_temp/{species}/{region}/{sample}_{hap}_{ref}_{region}_flank_genes_sort.bam",
        status = "regions/{species}/{sample}/{region}/{sample}_{hap}_{ref}_{region}_status.csv",
    output:
        bed = "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_flanking_genes.bed",
        bed_mod = temp("regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_coords.tmp"),
        bed_final = "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_region_coords.bed",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        # Convert to bed coords
        bedtools bamtobed -i {input.bam} > {output.bed}

        grep "closed" {input.status} |
        while IFS='\t' read -r flank_region flank_status genes contig; do
            if [[ "$flank_region" == "telomere" ]]; then
                # Can there be + and - strands together?
                strand=$(awk -v var="$contig" '$1==var {{print $6}}' {output.bed} | awk '!s[$0]++')
                echo $flank_region $flank_status $genes $contig $strand
                if [[ "$strand" == "+" ]]; then
                    end_coord=$(grep ":$contig\t" {input.sam} | grep "@" | awk -F'LN:' '{{print $2}}')
                    awk -v contig="$contig" -v endcoord="$end_coord" '$1==contig {{print $1"\t"$2"\t"endcoord"\t"$4"\t"$5"\t"$6}}' {output.bed} >> {output.bed_mod}
                else
                    awk '{{ print $1"\t0\t"$3"\t"$4"\t"$5"\t"$6}}' {output.bed} >> {output.bed_mod}
                fi
            else
            # Process closed flank region 
                # Smallest start and largest end
                awk -v var="$contig" '$1==var' {output.bed} |
                awk -v regin="{wildcards.region}" -F'\t' '{{
                  if (!start[$1] || $2 < start[$1]) {{
                    start[$1] = $2
                    region[$1] = regin
                    value[$1] = $5
                    strand[$1] = $6
                  }}
                  if (!end[$1] || $3 > end[$1]) {{
                    end[$1] = $3
                  }}
                }}
                END {{
                  for (key in start) {{
                    print key "\t" start[key] "\t" end[key] "\t" region[key] "\t" value[key] "\t" strand[key]
                  }}
                }}' >> {output.bed_mod}
            fi 
        done
        awk '!s[$2,$3]++' {output.bed_mod} > {output.bed_final}
        echo "Extract region with the flanking gene for {wildcards.sample} {wildcards.hap}"
        """

rule region_intact_minimap2:
    input:
        fa = "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_intact.fa",
        lib = lambda wildcards: REGIONS[wildcards.region]['library'],
    output:
        "map_temp/{species}/{region}/library/{sample}_{hap}_{ref}_{region}_intact.paf",
    params:
        minimap2 = lambda wildcards: REGIONS[wildcards.region]['minimap2'],
    threads: 6
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        # mapping lib to region
        minimap2 {params} -t {threads} {input.fa} {input.lib} > {output}
        echo "Mapping {wildcards.region} to {wildcards.sample} {wildcards.hap}"
        """

rule region_intact_cigar_info:
    input:
        "map_temp/{species}/{region}/library/{sample}_{hap}_{ref}_{region}_intact.paf",
    output:
        "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_gene_coords_paf.txt"
    shell:
        """
        python scripts/cigar_digest2.py {input} {output}
        echo "Processed paf file of {wildcards.sample} {wildcards.hap} for region {wildcards.region}"
        """

rule region_intact_contig:
    input:
        bed = "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_region_coords.bed",
        fa = "results/{sample}/scaffolds/{sample}_hybridPN_denovo_{hap}_{ref}.fa",
    output:
        fa = "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_intact.fa",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools getfasta -fi {input.fa} -bed {input.bed} -s | sed 's/:.*/_{wildcards.region}/g' > {output.fa}
        echo "Extract intact {wildcards.region} region with the flanking genes for {wildcards.sample} {wildcards.hap}"
        """

rule region_intact_gene_blast:
    input:
        coords = "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_gene_coords_paf.txt",
        fa = "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_intact.fa",
        lib = lambda wildcards: REGIONS[wildcards.region]['library']
    output:
        "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_gene_coords_paf.blast",
    threads: 10
    conda:
        "../envs/bbstools.yaml"
    params:
        blast = lambda wildcards: REGIONS[wildcards.region]['blast'],
    shell:
        """
        sed '/CIGAR/d' {input.coords} |
        awk '!s[$1]++ {{print $1}}' |
        for i in `cat`; do
            name=$(echo $i | sed 's/|.*//g')
            awk -v var="$i" '$1==var' {input.coords} | awk '{{print $2,$5,$6}}' | sed 's/ /\t/g' > regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name".bed"
            bedtools getfasta -fi {input.fa} -bed regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name".bed" > regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name"_region.fa" 
            seqkit grep -p "$i" {input.lib} > regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name"_lib.fasta"
            makeblastdb -in regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name"_lib.fasta" -dbtype nucl -logfile regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name"_lib.log"
            blastn -query regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name"_region.fa" -db regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name"_lib.fasta" {params.blast} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" >> {output}
            rm regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name".bed"  regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name"_region.fa"  regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name"_lib.fasta"*  regions/{wildcards.species}/{wildcards.sample}/{wildcards.region}/annotation/{wildcards.hap}/$name"_lib.log"
        done
        echo "Blast reference vs identified {wildcards.region} genes for {wildcards.sample} {wildcards.hap}"
        """

rule region_intact_join_blast_paf:
    input:
        info = "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_gene_coords_paf.txt",
        blast = "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_gene_coords_paf.blast",
    output:
        "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_analysis_paf.txt"
    shell:
        """
        python scripts/blast_extraction.py {input.info} {input.blast} -o {output}
        echo "Combined paf and blast files for {wildcards.region} region of {wildcards.sample} {wildcards.hap}"
        """

rule region_intact_process_blast_paf:
    input:
        "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_analysis_paf.txt"
    output:
        "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_query_groups.csv"
    params:
        "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_query_groups"
    shell:
        """
        python scripts/analyse_paf2.py {input} {params}
        echo "Process paf_blast files for {wildcards.region} region of {wildcards.sample} {wildcards.hap}"
        """

rule region_intact_filter_queries:
    input:
        "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_query_groups.csv"
    output:
        "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_query_groups_filtered.csv"
    params:
        "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_query_groups_filtered"
    shell:
        """
        python scripts/filter_query_ident.py {input} {params}
        echo "Filter for {wildcards.region} region of {wildcards.sample} {wildcards.hap}"
        """


rule region_intact_figure:
    input:
        "regions/{species}/{sample}/{region}/annotation/{hap}/{sample}_{hap}_{ref}_{region}_query_groups_filtered.csv"
    conda:
        "../envs/dna_viewer.yaml"
    output:
        "figures/{species}_{sample}_{hap}_{ref}_{region}.png"
    shell:
        """
        python scripts/DNA_viewer2_csv.py {input} {wildcards.species}_{wildcards.sample}_{wildcards.hap}_{wildcards.ref}_{wildcards.region}
        echo "Saving images for {wildcards.region} region of {wildcards.sample} {wildcards.hap}"
        """

rule closed_flanks:
    input:
        "figure/annotation_flow_{species}_{sample}_{hap}_{ref}_{region}.png"
    output:
        touch("regions/{species}/{sample}/{region}/closed/{sample}_{hap}_{ref}_{region}.done")
    shell:
        """
        echo "The intact {wildcards.region} region was analysed for sample {wildcards.sample} {wildcards.hap}"
        """

rule fragment_region:
    output:
        touch("regions/{species}/{sample}/{region}/fragmented/{sample}_{hap}_{ref}_{region}.done")
    shell:
        """
        echo "The region {wildcards.region} of sample {wildcards.sample} {wildcards.hap} was fragmented."
        """

rule missing_region:
    output:
        touch("regions/{species}/{sample}/{region}/missing/{sample}_{hap}_{ref}_{region}.done")

def check_gaps(wildcards):
    with checkpoints.check_flanking_genes.get(**wildcards).output[0].open() as f:
        for line in f:
            col1, col2, col3, col4 = line.strip().split("\t")
            if col2 == "closed":
                return "regions/{species}/{sample}/{region}/closed/{sample}_{hap}_{ref}_{region}.done"
            elif col2 == "fragmented":
                return "regions/{species}/{sample}/{region}/fragmented/{sample}_{hap}_{ref}_{region}.done"
            else:
                return "regions/{species}/{sample}/{region}/missing/{sample}_{hap}_{ref}_{region}.done"


rule check_region_annotation:
    input:
        check_gaps
    output:
        "regions/{species}/{sample}/{region}/flanking/{sample}_{hap}_{ref}_{region}_annotation.finish"
    shell:
        "cat {input} > {output}"

