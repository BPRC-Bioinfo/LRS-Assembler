# Filter based on blast results
# v0.0.8
# By Giang Le & Jaimy

import pandas as pd
import numpy as np

pd.set_option('display.max_rows', None)  # Show all rows
pd.set_option('display.max_columns', None)  # Show all columns (if necessary)
pd.set_option('display.width', None)  # Auto-adjust the width to display all columns
pd.set_option('display.max_colwidth', None)

def filter_percent_and_ref_len(df):
    # Percentage filter
    df_best_hits = df.sort_values(by=['percent', 'mismatch', 'roi_start'], ascending=[False, True, True])
    df_best_hits['vs_ref'] = abs(df_best_hits['align'] - df_best_hits['ref_len']) + df_best_hits['mismatch'] + df_best_hits.get('gap', 0) + df_best_hits.get('indel', 0)
    df_best_hits['align_percent'] = df_best_hits['align'] / df_best_hits['ref_len']
    # Filter hits below 70%
    df_best_hits = df_best_hits[df_best_hits['align_percent'] > 0.7]

    # Sort by highest blast_percent and lowest vs_ref, then group by roi_start and roi_end. Pick the first one.
    df_best_hits = (df_best_hits.sort_values(["percent", "vs_ref"], ascending=[False, True]).groupby(["roi_start", "roi_end"], as_index=False).first()) 
    
    return df_best_hits

def cluster_overlaps_per_strand(df_best_hits):
    """ Final result is a nested dictonary of strand and clusters"""
    clusters_by_strand = {}
    
    #Future groupby contigs ???
    for strand, group in df_best_hits.groupby('strand'):
        group_sorted = group.sort_values('roi_start').reset_index(drop=True)

        clusters = []   
        current_cluster = [] 
        current_cluster_end = None
        
        for idx, row in group_sorted.iterrows():
            hit = row.to_dict()

            if not current_cluster:
                current_cluster.append(hit)
                current_cluster_end = hit['roi_end']
            else:
                if hit['roi_start'] <= current_cluster_end:
                    current_cluster.append(hit)
                    current_cluster_end = max(current_cluster_end, hit['roi_end'])
                else:
                    clusters.append(current_cluster)
                    current_cluster = [hit]
                    current_cluster_end = hit['roi_end']
        
        if current_cluster:
            clusters.append(current_cluster)
        
        for cluster in clusters:
            cluster.sort(key=lambda hit: (hit['vs_ref'], -hit['percent']))
        
        cluster_dict = {i+1: cluster for i, cluster in enumerate(clusters)}
        clusters_by_strand[strand] = cluster_dict

    return clusters_by_strand

def group_overlapping_hits(hits):
    hits_sorted = sorted(hits, key=lambda h: h['roi_start'])
    groups = []
    current_group = []
    current_group_end = None
    for hit in hits_sorted:
        if not current_group:
            current_group = [hit]
            current_group_end = hit['roi_end']
        else:
            if hit['roi_start'] <= current_group_end:
                current_group.append(hit)
                current_group_end = max(current_group_end, hit['roi_end'])
            else:
                groups.append(current_group)
                current_group = [hit]
                current_group_end = hit['roi_end']
    if current_group:
        groups.append(current_group)
    
    return groups

def remove_contained(intervals):
    """Drop any interval that fully contains another."""
    return [
        I for I in intervals
        if not any(
            (I.roi_start <= J.roi_start and I.roi_end >= J.roi_end and I != J)
            for J in intervals
        )
    ]

def max_non_overlapping(intervals):
    """Greedily pick the max set of non‐overlapping intervals."""
    # sort by end
    sorted_iv = sorted(intervals, key=lambda x: x.roi_end)
    result = []
    last_end = -float('inf')
    for iv in sorted_iv:
        if iv.roi_start > last_end:
            result.append(iv)
            last_end = iv.roi_end
    return result

def cluster_hits_by_overlap(hits):
    """Group hits into clusters whenever their [roi_start, roi_end] intervals overlap."""
    # Union‑find setup
    n = len(hits)
    parent = list(range(n))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    # link any two overlapping hits
    for i in range(n):
        for j in range(i+1, n):
            hi, hj = hits[i], hits[j]
            if hi['roi_start'] <= hj['roi_end'] and hj['roi_start'] <= hi['roi_end']:
                union(i, j)

    # collect clusters by root parent
    clusters = {}
    for idx, hit in enumerate(hits):
        root = find(idx)
        clusters.setdefault(root, []).append(hit)
    return list(clusters.values())

def select_best_per_cluster(hits):
    """
    For each cluster of overlapping hits, pick the hit with:
      1) if tie, largest percent
      2) smallest vs_ref      
    """
    clusters = cluster_hits_by_overlap(hits)
    best_hits = []
    for cluster in clusters:
#        best = min(cluster, key=lambda h: (h['vs_ref'], -h['percent']))
        best = min(cluster, key=lambda h: (-h['percent'], h['vs_ref']))
        best_hits.append(best)
    return best_hits

def max_nonoverlapping(genes, start, end):
    """
    Given a list of gene-dicts with 'roi_start'/'roi_end',
    returns the maximal non-overlapping subset (greedy by end).
    """
    # 1) keep only those entirely inside [start, end]
    in_window = [g for g in genes
                 if g['roi_start'] >= start and g['roi_end'] <= end]
    # 2) sort by their end coordinate
    in_window.sort(key=lambda g: g['roi_end'])
    # 3) pick greedily
    selected = []
    last = start
    for g in in_window:
        if g['roi_start'] >= last:
            selected.append(g)
            last = g['roi_end']
    return selected

def process_hits_cDNA(hits):
    """
    - If any hits have vs_ref == 0:  do your zero‐vs_ref non‑overlap selection,
      then cluster the REMAINING hits and pick best-per-cluster.
    - If none have vs_ref == 0: cluster ALL hits and pick best-per-cluster.
    Returns a LIST of hit‑dicts.
    """
    if not hits:
        return []

    vs_ref_zero = [h for h in hits if h['vs_ref'] == 0]
#    vs_ref_zero = [h for h in hits if h['vs_ref'] == 0 and h['percent'] == 100.0]

    if vs_ref_zero:
        groups = group_overlapping_hits(vs_ref_zero)
#        print (groups)

        results = []
        for idx, group in enumerate(groups, start=1):
            if not group:
                print(f"Group {idx}: empty")
                continue

            roi_start = min(g['roi_start'] for g in group)
            roi_end   = max(g['roi_end']   for g in group)
            chosen = max_nonoverlapping(group, roi_start, roi_end)
            
            results.extend(chosen)           

        current_best = pd.DataFrame(results)

        # coords of current best
        intervals = current_best[['roi_start','roi_end']].to_records(index=False)

        filtered = remove_contained(intervals)
        selected = max_non_overlapping(filtered)

        selected_df = (
            current_best
            .set_index(['roi_start','roi_end'])
            .loc[[(iv.roi_start, iv.roi_end) for iv in selected]]
            .reset_index()
        )

        # remove anything overlapping those selected zeros
        sel_list = selected_df.to_dict('records')

        remaining = [
            h for h in hits
            if not any(
                h['roi_start'] <= s['roi_end'] and 
                h['roi_end']   >= s['roi_start']
                for s in sel_list
            )
        ]

        if remaining:
            remainingz = pd.DataFrame(remaining)
            remainingz = remainingz.sort_values("roi_start")

            # cluster & pick best from remaining
            best_from_clusters = select_best_per_cluster(remaining)

            # combine the two sets of winners
            combined = pd.concat([selected_df, pd.DataFrame(best_from_clusters)],
                                 ignore_index=True)
        else:
            combined =current_best.copy()
    else:
        best_from_clusters = select_best_per_cluster(hits)
        combined = pd.DataFrame(best_from_clusters)

    combined = combined.sort_values('roi_start').reset_index(drop=True)

    return combined.to_dict('records')

def process_hits_gDNA(hits):

    if not hits:
        return []

    result = []
    vs_ref_zero = [hit for hit in hits if hit['vs_ref'] == 0]
    
    if vs_ref_zero:
        groups = group_overlapping_hits(vs_ref_zero)
        current_best = [
            sorted(group, key=lambda h: (h['vs_ref'], -h['percent']))[0]
            for group in groups
        ]
        result.extend(current_best)
        
        remaining = []
        for hit in hits:
            if any(hit['roi_start'] <= best['roi_end'] and hit['roi_end'] >= best['roi_start']
                   for best in current_best):
                continue
            remaining.append(hit)
        
        # Recurse on the remaining hits.
        result.extend(process_hits_gDNA(remaining))
        return result
    else:
        best = sorted(hits, key=lambda h: (h['vs_ref'], -h['percent']))[0]
        result.append(best)
        remaining = [
            hit for hit in hits
            if not (hit['roi_start'] <= best['roi_end'] and hit['roi_end'] >= best['roi_start'])
        ]

        # Recurse on the remaining hits.
        result.extend(process_hits_gDNA(remaining))
        return result

# This function could be combined when pass lib type in 
# Check what is going on per example
        
def filtering_best_hits(clusters_by_strand, lib_type):
    """ Filter the best hit depending on lib input"""
    filtered_clusters = {}
    
    for strand, clusters in clusters_by_strand.items():
        filtered_clusters[strand] = {}
        for cluster_num, hits in clusters.items():
            if len(hits) == 1:
                best_hit = hits[0]
            else:
                if "gDNA" in lib_type:
                    best_hit = process_hits_gDNA(hits)
                elif "cDNA" in lib_type:
                    best_hit = process_hits_cDNA(hits)
            filtered_clusters[strand][cluster_num] = best_hit
   
    return filtered_clusters

def clusters_to_dataframe(filtered_clusters):
    records = []

    for strand, clusters in filtered_clusters.items():
        for cluster_num, best_hit in clusters.items():
            if isinstance(best_hit, list):
                for hit in best_hit:
                    record = hit.copy()
                    record['strand'] = strand
                    records.append(record)
            elif isinstance(best_hit, dict):
                record = best_hit.copy()
                record['strand'] = strand
                records.append(record)
            else:
                raise ValueError(f"Unexpected type {type(best_hit)} for cluster {cluster_num} in strand {strand}")

    return pd.DataFrame(records)

def assign_suffix(group):
    counter = 1
    suffixes = []

    for _, row in group.iterrows():
        if row['vs_ref'] == 0 and row['percent'] == 100:
            suffixes.append("")
        else:
            suffixes.append(f"_like{counter}")
            counter += 1
    group = group.copy() 
    group['gene_suffix'] = suffixes
    return group

def process_gene_names(df):
    df = df.copy()
    df['base_gene'] = df['gene_group'].str.split('|').str[0].str.replace(r'_group.*', '', regex=True)
    
    df_sorted = df.sort_values(by=['vs_ref', 'percent'], ascending=[True, False])
#    df_with_suffix = df_sorted.groupby('base_gene', group_keys=False).apply(assign_suffix)
    df_with_suffix = (
    df_sorted
      .groupby('base_gene', group_keys=True)
      .apply(assign_suffix, include_groups=False)
      .reset_index(level='base_gene')
    )
    df_with_suffix['gene_name'] = df_with_suffix['base_gene'] + df_with_suffix['gene_suffix'] + "_p" + df_with_suffix['percent'].astype(str)    
    df_with_suffix = df_with_suffix.drop(columns=['base_gene', 'gene_suffix', 'gene_group'])    
    df_with_suffix = df_with_suffix.sort_values('roi_start')

    return df_with_suffix

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Filter DataFrame based on aligned length and other criteria.")
    parser.add_argument("-i", "--input", help="Path to the input CSV file")
    parser.add_argument("-l", "--lib", help="Library type, cDNA or gDNA")
    parser.add_argument("-o", "--output", help="Path to the output Excel file")

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep = '\t') 

    filtered_df = filter_percent_and_ref_len(df)
    strand_cluster = cluster_overlaps_per_strand(filtered_df)
    filtered_cluster = filtering_best_hits(strand_cluster,args.lib)

    final_data = clusters_to_dataframe(filtered_cluster)
    final_data = final_data.sort_values('roi_start')
    final_data = process_gene_names(final_data)

    # Coords correction
    if "gDNA" in args.lib:
        final_data['roi_start'] = final_data['roi_start'] + 1
        final_data['ref_start'] = final_data['ref_start'] + 1
        final_data = final_data[['gene_name','roi_start','roi_end','percent','align','align_percent','mismatch','gap','ref_name','ref_start','ref_end','ref_len','vs_ref','contig','strand']].copy()
        mask = (
            final_data[['mismatch','gap','vs_ref']].ne(0).any(axis=1)
            |
            final_data['percent'].ne(100)
        )
    else:
        final_data = final_data[['gene_name','roi_start','roi_end','percent','align','align_percent','mismatch','indel','ref_name','ref_start','ref_end','ref_len','vs_ref','contig','strand']].copy()
        mask = (
            final_data[['mismatch','indel','vs_ref']].ne(0).any(axis=1)
            |
            final_data['percent'].ne(100)
        )
    final_data['ref_name'] = final_data['ref_name'].str.split('|').str[0]
    final_data['align_percent'] = final_data['align_percent'] * 100
    final_data['status'] = np.where(mask, 'novel', 'known')
    status = final_data.pop('status')

    loc = final_data.columns.get_loc('roi_end') + 1
    final_data.insert(loc, 'status', status)

#    print (final_data)
    final_data.to_csv(f"{args.output}.csv", sep='\t', float_format='%.2f' ,index=False)

    bed_file = final_data[['contig','roi_start','roi_end','gene_name','strand']].copy()
    bed_file['qual'] = 0
    cols = ['contig', 'roi_start', 'roi_end', 'gene_name', 'qual', 'strand']
    bed_file = bed_file[cols]
    bed_file.to_csv(f"{args.output}.bed", sep='\t' ,index=False, header=False)

