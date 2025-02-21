"""
Script combine blast results with paf file.

Author: Giang & Jaimy
Version: 0.0.4
Date: 2025-02-21
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
        chr_start = int(rest[3])
        chr_end = int(rest[4])

        return ref_name, int(ref_length), ref_start, ref_end, chr_start, chr_end, rest
    except ValueError:
        return None, None, None, None, None, None, None,

def group_by_gene(info_df):
    groups = {}
    
    for line in info_df:
        ref_name, ref_length, ref_start, ref_end, chr_start, chr_end, original_line = parse_line(line)

        if ref_name is not None:
            groups.setdefault(ref_name, []).append((ref_length, ref_start, ref_end, chr_start, chr_end, original_line))

    return groups

def individual_gene(genes):
    gene_groups = {}

    for ref_name, group_lines in genes.items():
        if ref_name is not None:
            current_ref_start = None
            subgroup = 0
            gene_groups[ref_name] = {}

            for i, (length, ref_start_coord, ref_end_coord, chr_start_coord, chr_end_coord, original_line) in enumerate(group_lines):
                if current_ref_start is None:
                    current_ref_start = ref_start_coord
                    current_ref_end = ref_end_coord
                    subgroup += 1
                    gene_groups[ref_name][subgroup] = [[original_line, 0 ,length]]
                elif ref_start_coord == current_ref_start or ref_start_coord < current_ref_start:
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

def join_blast_data(gene_groups, blast_df):

    result_data = []
    for ref_name, group_lines in gene_groups.items():
        for subgroup, rows in group_lines.items():
            for line in rows:
                contig = line[0][0]
                ref_start, ref_end, chr_start, chr_end, cigar, minimap_matches, minimap_inserts, minimap_deletions, strand = line[0][1:10]

                minimap_overlap = line[-2]
                blast_chr_coords = f"{contig}:{chr_start}-{chr_end}"
                seq_len = line[-1]
                ref_name_len = f"{ref_name}|{seq_len}"
                
                filtered_df = blast_df[(blast_df[0] == blast_chr_coords) & (blast_df[1] == ref_name_len)]
                
                blast_hits, blast_percent, blast_align, blast_overlap, blast_mismatch, blast_gap = [pd.NA] * 6
                
                if not filtered_df.empty:
                    first_evalue = filtered_df.iloc[0, 10]

                    if first_evalue != 0:
                        row = filtered_df.iloc[0]
                        blast_percent, blast_align, blast_mismatch, blast_gap, blast_ref_start, blast_ref_end, blast_chr_start, blast_chr_end = row[2:10]
                        overlap_sum = 0
                    else:
                        zero_evalue_rows = filtered_df[filtered_df[10] == 0].sort_values(by=6)
                        smallest_start = zero_evalue_rows[6].min()
                        largest_end = zero_evalue_rows[7].max()

                        # Apply the filter to remove rows that fall strictly within the largest fragment
                        within_coords_filter = zero_evalue_rows[~zero_evalue_rows.apply(falls_within_largest, axis=1, args=(smallest_start, largest_end))]

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
                                    difference = blast_ref_end - current_blast_ref_start
                                    blast_overlap.append(difference)
                            blast_ref_start = current_blast_ref_start
                            blast_ref_end = current_blast_ref_end

                        overlap_sum = 0 if not blast_overlap else sum(blast_overlap)
                        blast_align = (within_coords_filter[7] - within_coords_filter[6] - 1).sum()
                        blast_mismatch, blast_gap = within_coords_filter[4].sum(), within_coords_filter[5].sum()

                    blast_hits = len(filtered_df)
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
                    "chr_start": int(chr_start),
                    "chr_end": int(chr_end),
                    "minimap_matches": int(minimap_matches),
                    "minimap_inserts": int(minimap_inserts),
                    "minimap_deletions": int(minimap_deletions),
                    "minimap_overlap": int(minimap_overlap),
                    "strand": strand,
                    "blast_hits": blast_hits,
                    "blast_percent": blast_percent,
                    "blast_align": blast_align,
                    "blast_overlap": blast_overlap,
                    "blast_mismatch": blast_mismatch,
                    "blast_gap": blast_gap,
                    "ref_len": int(seq_len)
                })

    paf_blast_df = pd.DataFrame(result_data)
    paf_blast_df = paf_blast_df.sort_values('chr_start')

    paf_blast_regroup_df = paf_blast_df.sort_values(['ref_name', 'gene_group', 'ref_start']).reset_index(drop=True)
    new_group = paf_blast_regroup_df['gene_group'].max() + 1

    prev = None
    # Regroup if fragments larger than 70%
    # Maybe do condition if total fragments is larger than 150%
    for idx, row in paf_blast_regroup_df.iterrows():
        if prev is not None:
            if row['ref_name'] == prev['ref_name'] and row['gene_group'] == prev['gene_group']:
                if row['chr_start'] != prev['chr_end']:
                    blast_align_val = row['blast_align'] if pd.notna(row['blast_align']) else 0
                    # Threshold for splitting gene group is 0.7
                    if (blast_align_val / row['ref_len']) > 0.7:
                        paf_blast_regroup_df.at[idx, 'gene_group'] = new_group
                        paf_blast_regroup_df.at[idx, 'minimap_overlap'] = 0
                        new_group += 1
        prev = row

    return paf_blast_regroup_df

def process_individual_gene(paf_blast_df):
    final_data = []
    for (ref_name, gene_group), group in paf_blast_df.groupby(['ref_name', 'gene_group']):
        sorted_group = group.sort_values(by='ref_start')

        ref_start_fin = sorted_group['ref_start'].iloc[0]
        ref_end_fin = sorted_group['ref_end'].iloc[-1]
        chr_start_fin = sorted_group['chr_start'].iloc[0]
        chr_end_fin = sorted_group['chr_end'].iloc[-1]
        strand = sorted_group['strand'].iloc[0]
        ref_len = sorted_group['ref_len'].iloc[0]

        minimap_matches_fin = sorted_group['minimap_matches'].sum()
        minimap_inserts_fin = sorted_group['minimap_inserts'].sum()        
        minimap_deletions_fin = sorted_group['minimap_deletions'].sum()
        minimap_overlap_fin = sorted_group['minimap_overlap'].sum()

        num_nan_rows = sorted_group['blast_percent'].isna().sum()
        blast_percent_fin = sorted_group['blast_percent'].dropna().mean()
        blast_align_fin = sorted_group['blast_align'].dropna().sum() - sorted_group['minimap_overlap'].dropna().sum() - sorted_group['blast_overlap'].dropna().sum()
        blast_align_percent_fin = blast_align_fin / ref_len * 100
        blast_mismatch_fin = sorted_group['blast_mismatch'].dropna().sum()
        blast_gap_fin = sorted_group['blast_gap'].dropna().sum()

        final_data.append({
            "gene_group": f"{ref_name.split('|')[0]}_group{gene_group}|{ref_len}",
            "ref_name": ref_name,
            "ref_start": int(ref_start_fin),
            "ref_end": int(ref_end_fin),
            "chr_start": int(chr_start_fin),
            "chr_end": int(chr_end_fin),
            "minimap_matches": int(minimap_matches_fin),
            "minimap_inserts": int(minimap_inserts_fin),
            "minimap_deletions": int(minimap_deletions_fin),
            "minimap_overlap": int(minimap_overlap_fin),
            "strand": strand,
            "blast_percent": f"{blast_percent_fin:.3f}",
            "blast_align": f"{blast_align_percent_fin:.3f}",
            "blast_mismatch": int(blast_mismatch_fin),
            "blast_gap": int(blast_gap_fin),
            "ref_len": int(ref_len)
        })

    paf_blast_final = pd.DataFrame(final_data)
    paf_blast_final = paf_blast_final.sort_values('chr_start')
    return paf_blast_final

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process info and blast files.")
    parser.add_argument("-i", "--info", help="Path to the info file.")
    parser.add_argument("-b", "--blast", help="Path to the blast file.")
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
    paf_blast_df = join_blast_data(multiple_genes, blast_df)
    genes_paf_blast = process_individual_gene(paf_blast_df)

#    print (genes_paf_blast)
    genes_paf_blast.to_csv(args.output, sep="\t", index=False)


'''
#    print (paf_blast_regroup_df[paf_blast_regroup_df['ref_name'].str.contains("LILRB2-like_4-1-1")])
#    print (genes_paf_blast[genes_paf_blast['ref_name'].str.contains("LILRB2-like_2-1-1")])
'''
