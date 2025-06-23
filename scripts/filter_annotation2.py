# Filter based on blast results
# v0.0.9
# By Giang Le & Jaimy

import pandas as pd
import numpy as np

pd.set_option('display.max_rows', None)  # Show all rows
pd.set_option('display.max_columns', None)  # Show all columns (if necessary)
pd.set_option('display.width', None)  # Auto-adjust the width to display all columns
pd.set_option('display.max_colwidth', None)

def filter_percent_and_ref_len(df, lib_type):
    # Percentage filter
    if "protein" not in lib_type:
        df_best_hits = df.sort_values(by=['percent', 'mismatch', 'roi_start'], ascending=[False, True, True])
    else:
        df_best_hits = df.sort_values(by=['percent', 'roi_start'], ascending=[False, True])

    df_best_hits['vs_ref'] = abs(df_best_hits['align'] - df_best_hits['ref_len']) + df_best_hits.get('mismatch',0) + df_best_hits.get('gap', 0) + df_best_hits.get('indel', 0)
    df_best_hits['align_percent'] = df_best_hits['align'] / df_best_hits['ref_len']

    # Filter hits below 70%
    df_best_hits = df_best_hits[df_best_hits['align_percent'] > 0.7]

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
#    print (sorted_iv)
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
    For each cluster of overlapping hits, pick:
      - among hits with percent == 100.0:
          • if they all share the same align → take them all
          • otherwise → take the one with the highest align
      - otherwise (no percent==100): the hit with highest percent (tie-break by smallest vs_ref)
    """
    best_hits = []
    clusters = cluster_hits_by_overlap(hits)

    for cluster in clusters:
        # normalize types
        for h in cluster:
            h['percent'] = float(h['percent'])
            h['vs_ref']  = int(h['vs_ref'])
            h['align_percent']   = float(h['align_percent'])

        # look for perfect matches
        perfect = [h for h in cluster if h['percent'] == 100.0]
#        print (perfect)
        if perfect:
            # check if all their align values are identical
            align_vals = {h['align'] for h in perfect}
#            print (align_vals)
            if len(align_vals) == 1:
                # all ties → keep them all
                best_hits.extend(perfect)
            else:
                # pick the one perfect hit with the highest align
                # find the highest align_percent
                max_ap = max(h['align_percent'] for h in perfect)
                # pick *all* perfect hits with that align_percent
                winners = [h for h in perfect if h['align_percent'] == max_ap]
                best_hits.extend(winners)

        else:
            # no 100% hits → fall back to highest percent, then smallest vs_ref
            max_pct = max(h['percent'] for h in cluster)
#            print (max_pct)
            candidates = [h for h in cluster if h['percent'] == max_pct]
            best = min(candidates, key=lambda h: h['vs_ref'])
            best_hits.append(best)

    return best_hits


def max_nonoverlapping(genes, start, end):
    """
    Given a list of gene-dicts with 'roi_start'/'roi_end',
    returns the maximal non-overlapping subset (greedy by end).
    If two candidates have the same roi_end *and* the same vs_ref/percent,
    the shorter interval (smaller roi_end–roi_start) will be picked first.
    """
    # 1) keep only those entirely inside [start, end]
    in_window = [g for g in genes
                 if g['roi_start'] >= start and g['roi_end'] <= end]

    # 2) sort by:
    #    (a) end coordinate
    #    (b) vs_ref      (lower is “better”)
    #    (c) percent     (higher is “better”, so we use -percent)
    #    (d) size        (smaller is “better”)
    in_window.sort(key=lambda g: (
        g['roi_end'],
        g.get('vs_ref', 0),
        -g.get('percent', 0.0),
        (g['roi_end'] - g['roi_start'])
    ))
#    print (in_window)
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
    Recursively select hits to fill in the whole ROI:
    1) If any hits have vs_ref == 0, cluster those and pick non-overlapping bests.
       Else cluster ALL hits and pick best-per-cluster.
    2) Remove selected hits from the pool.
    3) Recurse on the remainder until empty.
    Returns a LIST of hit-dicts, sorted by roi_start.
    """
    if not hits:
        return []
#    print (hits)
    # 1) choose this round’s winners
    vs_ref_zero = [h for h in hits if h['vs_ref'] == 0]

#    vs_ref_zero_df = pd.DataFrame(vs_ref_zero)
#    print (vs_ref_zero_df)

    if vs_ref_zero:
        # cluster zeros and pick non-overlapping bests
        groups = group_overlapping_hits(vs_ref_zero)

        interim = []
        for group in groups:
            if not group:
                continue
            start = min(h['roi_start'] for h in group)
            end   = max(h['roi_end']   for h in group)
            interim.extend(max_nonoverlapping(group, start, end))
        # filter contained intervals and pick global max‐nonoverlapping
        interim_df = pd.DataFrame(interim)

        ivs = interim_df[['roi_start','roi_end']].to_records(index=False)
#        print (interim_df['contig'].unique(), len(ivs))
        filtered_coords = remove_contained(ivs)
#        print (interim_df['contig'].unique(), len(filtered_coords))
        best_ivs = max_non_overlapping(filtered_coords)
#        print (interim_df['contig'].unique(), len(best_ivs))
        winners = (
            interim_df
            .set_index(['roi_start','roi_end'])
            .loc[[(iv.roi_start, iv.roi_end) for iv in best_ivs]]
            .reset_index()
            .to_dict('records')
        )
    else:
        # no zeros: just pick best per cluster over all hits
        winners = select_best_per_cluster(hits)
#    print ("(*_*)" * 10)
    
    winners_df = pd.DataFrame(winners)
#    print (winners_df)
#    print ("(T_T)" * 10)
#    print (hits)
    # 2) subtract overlapping hits
    remaining = [
        h for h in hits
        if not any(
            h['roi_start'] <= w['roi_end'] and
            h['roi_end']   >= w['roi_start']
            for w in winners
        )
    ]

#    print (remaining)
    # 3) recurse on what's left
    return sorted(
        winners + process_hits_cDNA(remaining),
        key=lambda h: h['roi_start']
    )

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
        # 1) pick a representative “best” to read off vs_ref, percent, align
        best = sorted(
            hits,
            key=lambda h: (h['vs_ref'], -h['percent'], h['align_percent'])
        )[0]
        target_vs_ref = best['vs_ref']
        target_pct    = best['percent']
        target_align  = best['align_percent']

        # 2) collect all hits that tie on those three fields
        best_hits = [
            h for h in hits
            if (h['vs_ref']   == target_vs_ref
                and h['percent'] == target_pct
                and h['align_percent']   == target_align)
        ]

        # add them all to the result
        result.extend(best_hits)

        # 3) filter out any hit that overlaps _any_ of the chosen best_hits
        remaining = [
            h for h in hits
            if not any(
                h['roi_start'] <= bh['roi_end'] and
                h['roi_end']   >= bh['roi_start']
                for bh in best_hits
            )
        ]

        # 4) recurse on whatever’s left
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
#            print (hits)
            if len(hits) == 1:
                best_hit = hits[0]
            else:
                if "cDNA" in lib_type:
                    best_hit = process_hits_cDNA(hits)
                else:
                    best_hit = process_hits_gDNA(hits)
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
    df_with_suffix = (
    df_sorted
      .groupby('base_gene', group_keys=True)
      .apply(assign_suffix, include_groups=False)
      .reset_index(level='base_gene')
    )
    df_with_suffix['gene_name'] = (
    df_with_suffix['base_gene']
    + df_with_suffix['gene_suffix']
    + "_p"
    + df_with_suffix['percent'].map("{:.2f}".format)
    )   
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
#    print (df)

    filtered_df = filter_percent_and_ref_len(df,args.lib)
#    print (filtered_df)
    strand_cluster = cluster_overlaps_per_strand(filtered_df)

#    print (strand_cluster)
    
    filtered_cluster = filtering_best_hits(strand_cluster, args.lib)
#    print (filtered_cluster)

    final_data = clusters_to_dataframe(filtered_cluster)
    print (final_data)

    matching_cols = ['roi_start', 'roi_end', 'percent', 'strand', 'vs_ref', 'ref_len', 'contig']

    merged = filtered_df.merge(
        final_data[matching_cols + ['ref_name']],  # include ref_name to exclude it later
        on=matching_cols,
        how='inner',
        suffixes=('', '_final')
    )

    # Step 2: Remove rows that have the same ref_name (only keep other hits)
    alternative_hits = merged[merged['ref_name'] != merged['ref_name_final']].drop(columns=['ref_name_final'])

    alternative_hits = alternative_hits[final_data.columns]
    duplicate_keys = final_data[['ref_name', 'roi_start', 'roi_end']].drop_duplicates()
    alternative_hits = alternative_hits.merge(
        duplicate_keys,
        on=['ref_name', 'roi_start', 'roi_end'],
        how='left',
        indicator=True
    ).query('_merge == "left_only"').drop(columns=['_merge'])

    # Step 3: Align and concat
    alternative_hits = alternative_hits[final_data.columns]
    final_data = pd.concat([final_data, alternative_hits], ignore_index=True).drop_duplicates()
    final_data = final_data.sort_values(['roi_start', 'roi_end']).reset_index(drop=True)
    final_data = process_gene_names(final_data)

    # Coords correction
    
    if "cDNA" in args.lib:
        final_data['roi_start'] = final_data['roi_start'] - 1
        final_data['ref_start'] = final_data['ref_start'] - 1

        final_data = final_data[['gene_name','roi_start','roi_end','percent','align','align_percent','mismatch','indel','ref_name','ref_start','ref_end','ref_len','vs_ref','contig','strand']].copy()
        mask = (
            final_data[['mismatch','indel','vs_ref']].ne(0).any(axis=1)
            |
            final_data['percent'].ne(100)
        )
    elif "gDNA" in args.lib:
        final_data['roi_start'] = final_data['roi_start'] + 1
        final_data['ref_start'] = final_data['ref_start'] + 1

        final_data = final_data[['gene_name','roi_start','roi_end','percent','align','align_percent','mismatch','gap','ref_name','ref_start','ref_end','ref_len','vs_ref','contig','strand']].copy()
        mask = (
            final_data[['mismatch','gap','vs_ref']].ne(0).any(axis=1)
            |
            final_data['percent'].ne(100)
        )
    else:
        final_data = final_data[['gene_name','roi_start','roi_end','percent','align','align_percent','ref_name','ref_start','ref_end','ref_len','vs_ref','contig','strand']].copy()
        mask = (
            final_data['vs_ref'].ne(0)
            |
            final_data['percent'].ne(100)
        )
        
    final_data['ref_name'] = final_data['ref_name'].str.split('|').str[0]
    final_data['align_percent'] = final_data['align_percent'] * 100
    final_data['status'] = np.where(mask, 'novel', 'known')
    status = final_data.pop('status')
    final_data['library'] = args.lib

    loc = final_data.columns.get_loc('roi_end') + 1
    final_data.insert(loc, 'status', status)

    print (final_data)
    final_data.to_csv(f"{args.output}.csv", sep='\t', float_format='%.2f' ,index=False)

    bed_file = final_data[['contig','roi_start','roi_end','gene_name','strand']].copy()
    bed_file['qual'] = 0
    cols = ['contig', 'roi_start', 'roi_end', 'gene_name', 'qual', 'strand']
    bed_file = bed_file[cols]
    bed_file.to_csv(f"{args.output}.bed", sep='\t' ,index=False, header=False)

