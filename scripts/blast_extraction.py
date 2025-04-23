"""
Script combine blast results with paf file.

Author: Giang & Jaimy
Version: 0.0.7
"""

import pandas as pd
import argparse

pd.set_option('display.max_columns', None) 

def parse_line(line):
    try:
        ref_info, *rest = line.strip().split("\t") 
        ref_name, ref_length = ref_info.split("|") 
        ref_start = int(rest[1])
        ref_end = int(rest[2])
        roi_start = int(rest[3])
        roi_end = int(rest[4])
        strand = rest[9]

        return ref_name, int(ref_length), ref_start, ref_end, roi_start, roi_end, strand, rest
    except ValueError:
        return None, None, None, None, None, None, None, None

def group_by_gene(info_df):
    genes = {}
    
    for line in info_df:
        ref_name, ref_length, ref_start, ref_end, roi_start, roi_end, strand, original_line = parse_line(line)

        if ref_name is not None:
            genes.setdefault(ref_name, []).append((ref_length, ref_start, ref_end, roi_start, roi_end, strand, original_line))

    return genes

def individual_gene(genes):

    for gene, entries in genes.items():
        genes[gene] = sorted(entries, key=lambda x: x[1])

    gene_groups = {}

    for ref_name, group_lines in genes.items():
        if ref_name is not None:
            gene_groups[ref_name] = {}
            subgroup = 0  
            
            for strand in ['+', '-']:
                strand_lines = [line for line in group_lines if line[5] == strand]
                if not strand_lines:
                    continue 

                current_ref_start = None
                for (length, ref_start_coord, ref_end_coord, roi_start_coord, roi_end_coord, current_strand, original_line) in strand_lines:
                    if current_ref_start is None:
                        current_ref_start = ref_start_coord
                        current_ref_end = ref_end_coord
                        subgroup += 1
                        gene_groups[ref_name][subgroup] = [[original_line, 0, length]]
                    elif ref_start_coord == current_ref_start or ref_start_coord < current_ref_start:
                        current_ref_start = ref_start_coord
                        current_ref_end = ref_end_coord
                        subgroup += 1
                        gene_groups[ref_name][subgroup] = [[original_line, 0, length]]
                    else:
                        if current_ref_end > ref_start_coord:
                            overlap = current_ref_end - ref_start_coord
                            gene_groups[ref_name][subgroup].append([original_line, overlap, length])
                        else:
                            gene_groups[ref_name][subgroup].append([original_line, 0, length])
                        current_ref_start = ref_start_coord

    return gene_groups

def falls_within_largest(row, smallest_start, largest_end):
    start_in_range = smallest_start < row[6] < largest_end
    end_in_range = smallest_start < row[7] < largest_end
    return start_in_range and end_in_range

def merge_blast_data(gene_groups, blast_df):

    result_data = []
    for ref_name, group_lines in gene_groups.items():
        for subgroup, rows in group_lines.items():
            for line in rows:
                contig = line[0][0]
                ref_start, ref_end, roi_start, roi_end, cigar, minimap_matches, minimap_inserts, minimap_deletions, strand = line[0][1:10]

                minimap_overlap = line[-2]
                blast_chr_coords = f"{contig}:{roi_start}-{roi_end}"
                seq_len = line[-1]
                ref_name_len = f"{ref_name}|{seq_len}"
                
                filtered_df = blast_df[(blast_df[0] == blast_chr_coords) & (blast_df[1] == ref_name_len)]
                
                blast_hits, blast_percent, blast_align, blast_overlap, blast_mismatch, blast_gap = [pd.NA] * 6
                
                if not filtered_df.empty:
                    first_evalue = filtered_df.iloc[0, 10]

                    if first_evalue != 0:
                        row = filtered_df.iloc[0]
                        blast_percent, blast_align, blast_mismatch, blast_gap, blast_ref_start, blast_ref_end, blast_roi_start, blast_roi_end = row[2:10]
                        overlap_sum = 0
                    else:
                        zero_evalue_rows = filtered_df[filtered_df[10] == 0].sort_values(by=6) # sort by ref start
                        recorded_ranges = []
                        keep_indices = []

                        for idx, row in zero_evalue_rows.iterrows():
                            # Get the values from column 7 and 8 (1-indexed) i.e. columns 6 and 7 (0-indexed)
                            start_val = row.iloc[6]
                            end_val = row.iloc[7]
                            
                            current_range = (start_val, end_val)
                            remove = False
                            # Check if current_range is within any previously recorded range
                            for rec in recorded_ranges:
                                if current_range[0] >= rec[0] and current_range[1] <= rec[1]:
                                    remove = True
                                    break
                            # If not contained in any previous range, record it and keep the row
                            if not remove:
                                recorded_ranges.append(current_range)
                                keep_indices.append(idx)

                        # Create a new DataFrame with only the rows we want to keep
                        within_coords_filter = zero_evalue_rows.loc[keep_indices].reset_index(drop=True)
                        
                        blast_percent = within_coords_filter[2].mean()
                        blast_ref_start, blast_ref_end = None, None
                        blast_overlap = []

                        for index, blast_row in within_coords_filter.iterrows():
                            current_blast_ref_start = int(blast_row[6])
                            current_blast_ref_end = int(blast_row[7])

                            if blast_ref_start is None:
                                blast_ref_start = current_blast_ref_start
                                blast_ref_end = current_blast_ref_end
                            else:
                                if current_blast_ref_start < blast_ref_end:
                                    difference = blast_ref_end - (current_blast_ref_start - 1)
                                    blast_overlap.append(difference)
                            blast_ref_start = current_blast_ref_start
                            blast_ref_end = current_blast_ref_end

                        overlap_sum = 0 if not blast_overlap else sum(blast_overlap) 
                        blast_align = (within_coords_filter[7] - (within_coords_filter[6] - 1)).sum() - overlap_sum 
                        blast_mismatch, blast_gap = within_coords_filter[4].sum(), within_coords_filter[5].sum()

                    blast_percent = float(blast_percent)
                    blast_align = int(blast_align)
                    blast_overlap = int(overlap_sum)
                    blast_mismatch = int(blast_mismatch)
                    blast_gap = int(blast_gap)

                result_data.append({
                    "ref_name": ref_name_len,
                    "gene_group": subgroup,
                    "ref_start": int(ref_start),
                    "ref_end": int(ref_end),
                    "roi_start": int(roi_start),
                    "roi_end": int(roi_end),
                    "minimap_matches": int(minimap_matches),
                    "minimap_inserts": int(minimap_inserts),
                    "minimap_deletions": int(minimap_deletions),
                    "minimap_overlap": int(minimap_overlap),
                    "strand": strand,
                    "percent": blast_percent,
                    "align": blast_align,
                    "mismatch": blast_mismatch,
                    "gap": blast_gap,
                    "contig": contig,
                    "ref_len": int(seq_len)
                })

    paf_blast_df = pd.DataFrame(result_data)
    paf_blast_df = paf_blast_df.sort_values('roi_start')

    return paf_blast_df

def regroup_genes(paf_blast_df):
    clusters_list = []
    sum_cols = ['minimap_matches', 'minimap_inserts', 'minimap_deletions',
                'minimap_overlap', 'align', 'mismatch', 'gap']
    
    for (ref_name, gene_group), group in paf_blast_df.groupby(['ref_name', 'gene_group']):
        current_cluster = []
        
        if len(group) == 1:
            for idx, row in group.iterrows():
                current_cluster.append(row.copy())
            clusters_list.append(current_cluster)
        else:
            strand = group.iloc[0]['strand']
            if strand == '-':
                merge_condition = lambda last, curr: last['ref_start'] > curr['ref_end']
                boundary_field = 'ref_start'
            else:
                merge_condition = lambda last, curr: last['ref_end'] < curr['ref_start']
                boundary_field = 'ref_end'
            
            for idx, row in group.iterrows():
                row = row.copy()
                if not current_cluster:
                    current_cluster.append(row)
                else:
                    last_row = current_cluster[-1]
                    if merge_condition(last_row, row):
                        last_row[boundary_field] = row[boundary_field]
                        last_row['roi_end'] = row['roi_end']
                        for col in sum_cols:
                            last_row[col] += row[col]
                        last_row['percent'] = (last_row['percent'] + row['percent']) / 2
                        current_cluster[-1] = last_row
                    else:
                        clusters_list.append(current_cluster)
                        current_cluster = [row]
            if current_cluster:
                clusters_list.append(current_cluster)
    
    merged_rows = [cluster[-1] for cluster in clusters_list]
    regroup_genes = pd.DataFrame(merged_rows)

    return regroup_genes

def rename_group(regrouped_genes):
    regrouped_genes = regrouped_genes.drop('minimap_overlap', axis = 1)
    updated_df = regrouped_genes.copy()
    for ref_name, ref_group in regrouped_genes.groupby("ref_name"):

        current_max = ref_group["gene_group"].max()
        
        for gene_group, group in ref_group.groupby("gene_group"):
            if len(group) > 1:
                for i, idx in enumerate(group.index):
                    if i > 0:
                        current_max += 1
                        updated_df.at[idx, "gene_group"] = current_max
    
    updated_df["gene_group"] = updated_df.apply(
        lambda row: f"{row['ref_name'].split('|')[0]}_group{row['gene_group']}|{row['ref_len']}",
        axis=1
    )

    updated_df = updated_df.sort_values('roi_start')

    return updated_df



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process info and blast files.")
    parser.add_argument("-i", "--info", help="Path to the info file.")
    parser.add_argument("-b", "--blast", help="Path to the blast file.")
    parser.add_argument("-l", "--lib", help="Library either cDNA or gDNA")
    parser.add_argument("-o", "--output", default="output.txt", help="Path to the output file.")
    args = parser.parse_args()

    try:
        with open(args.info, 'r') as f:
            lines = f.readlines()
        blast_df = pd.read_csv(args.blast, sep="\t", header=None)
    except FileNotFoundError: # Redundant check after os.path.exists, but kept for robustness
        print("Error: One or both of the files specified were not found (even after initial existence check).")
        exit(1)
    except pd.errors.ParserError:
        print("Error: Could not parse one or both of the input files. Ensure they are valid TSV files.")
        exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading the input files: {e}")
        exit(1)

    genes = group_by_gene(lines)
    multiple_genes = individual_gene(genes)
    paf_blast_df = merge_blast_data(multiple_genes, blast_df)
    similar_join = regroup_genes(paf_blast_df)
    final_genes_paf_blast = rename_group(similar_join)
#    print (final_genes_paf_blast.to_string())
    final_genes_paf_blast.to_csv(args.output, sep="\t", index=False)


