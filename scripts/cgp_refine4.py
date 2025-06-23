import argparse
import pandas as pd
import sys
import os
import re
from os.path import commonprefix
from itertools import takewhile
from typing import List, Tuple


def parse_args():
    parser = argparse.ArgumentParser(
        description='Load, cluster, and analyze cDNA, gDNA and/or protein TSV files.'
    )
    parser.add_argument('--input', '-i', help='Path to the combined TSV file (or "None").')
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Prefix for output files.'
    )
    return parser.parse_args()

def load_tsv(path):
    """
    Helper to read a TSV into a DataFrame. Exits on failure.
    """
    try:
        return pd.read_csv(path, sep='\t', header=0, dtype=str)
    except Exception as e:
        sys.exit(f'Error: Could not read {path!r} as a TSV. Exception:\n  {e}')

def infer_output_prefix(input_paths):
    """
    If the user did not specify --output, use the common prefix of all
    basenames (without extension). E.g. ['A_cdna.tsv','A_gdna.tsv'] → 'A'.
    """
    basenames = [os.path.splitext(os.path.basename(p))[0] for p in input_paths]
    if not basenames:
        return None

    cp = commonprefix(basenames)
    cp = cp.rstrip('_-')
    return cp if cp else None

def main():
    args = parse_args()

    if not args.input or args.input == "None":
        sys.exit("Error: You must supply --input with a valid combined file.")

    if not os.path.isfile(args.input):
        sys.exit(f"Error: File not found: {args.input!r}")

    df = load_tsv(args.input)

    if 'library' not in df.columns:
        sys.exit("Error: Input file must contain a 'library' column.")

    library_kinds = set(df['library'].dropna().unique())

    if not library_kinds:
        sys.exit("Error: No valid 'library' types found in the input file.")

    out_prefix = args.output if args.output and args.output != "None" else infer_output_prefix([args.input]) or "combined"

    print("Detected the following library types:")
    for kind in sorted(library_kinds):
        print(f"  - {kind.upper()}")
    print(f"Using output prefix: {out_prefix!r}\n")

    # Dispatch to appropriate function
    if {'cDNA', 'gDNA', 'protein'} <= library_kinds:
        run_all_three(df[df['library'] == 'cDNA'],
                      df[df['library'] == 'gDNA'],
                      df[df['library'] == 'protein'],
                      out_prefix)
    elif {'cDNA', 'gDNA'} <= library_kinds:
        run_cdna_gdna(df[df['library'] == 'cDNA'],
                      df[df['library'] == 'gDNA'],
                      out_prefix)
    elif {'cDNA', 'protein'} <= library_kinds:
        run_cdna_protein(df[df['library'] == 'cDNA'],
                         df[df['library'] == 'protein'],
                         out_prefix)
    elif {'gDNA', 'protein'} <= library_kinds:
        run_gdna_protein(df[df['library'] == 'gDNA'],
                         df[df['library'] == 'protein'],
                         out_prefix)
    elif 'cDNA' in library_kinds:
        run_single_cdna(df[df['library'] == 'cDNA'], out_prefix)
    elif 'gDNA' in library_kinds:
        run_single_gdna(df[df['library'] == 'gDNA'], out_prefix)
    elif 'protein' in library_kinds:
        run_single_protein(df[df['library'] == 'protein'], out_prefix)
    else:
        sys.exit("Error: Unhandled combination of 'library' types.")

def cluster_intervals2(df, start_col, end_col, group_cols):
    df_copy = df.copy()
    df_copy['cluster_id'] = None
    for grp_vals, grp_df in df_copy.groupby(group_cols):
        sorted_idx = grp_df[start_col].sort_values().index
        clusters, current, current_end = [], [], None
        for idx in sorted_idx:
            s, e = grp_df.at[idx, start_col], grp_df.at[idx, end_col]
            if current_end is None or s > current_end:
                if current:
                    clusters.append(current)
                current, current_end = [idx], e
            else:
                current.append(idx)
                current_end = max(current_end, e)
        if current:
            clusters.append(current)
        for cid, cluster in enumerate(clusters, start=1):
            label = f"cluster_{'_'.join(map(str, grp_vals))}_{cid}"
            for idx in cluster:
                df_copy.at[idx, 'cluster_id'] = label
    return df_copy

# === UTILITIES ===
def lcp(a: str, b: str) -> str:
    """Longest common prefix between two strings."""
    return ''.join(x for x, _ in takewhile(lambda x: x[0] == x[1], zip(a, b)))

def prefix_similarity(a: str, b: str) -> float:
    """Return ratio of longest common prefix length over max string length."""
    prefix = lcp(a, b)
    return len(prefix) / max(len(a), len(b))

# === CLUSTER PROCESSING ===
def process_cluster2(cluster_df: pd.DataFrame):
    # BEST OVERALL HIT
    top_hit = _best_hit(cluster_df)
    best_name, best_src = top_hit['ref_name'], top_hit['source']
    sources = set(cluster_df['source'])

    opp = None
    # SINGLE-SOURCE CASE
    if len(sources) == 1:
        opp = _best_hit(cluster_df, exclude_source=best_src)
    
    # MULTI-SOURCE CASE
    opp = _best_hit(cluster_df, exclude_source=best_src)

    # PERFECT KNOWN HITS
    if set(top_hit['source']) == {'cDNA', 'gDNA', 'protein'} and \
       (top_hit['vs_ref'].astype(float) == 0).all() and \
       (top_hit['percent'].astype(float) == 100).all():
        names_list = top_hit['ref_name'].tolist()
        # All have same ref_name
        if len(set(names_list)) == 1:
            gene_name = names_list[0]
        else:
            # Pick gDNA if all 3 hits are identical match
            # gDNA most likely for changes
            gene_name = top_hit[top_hit['source']=='gDNA'].iloc[0]['ref_name']
        avg_pct, start, end = _cluster_metrics(top_hit, None)
        return gene_name, avg_pct, start, end, 'c+g+p'
    
    # Multiple top hits
    if len(top_hit) > 1:
        print (top_hit)
        # INTERNAL PREFIX SIMILARITY MATRIX
        internal_rows = []
        top_names = top_hit['ref_name'].tolist()

        for i, hi in top_hit.iterrows():
            for j, hj in top_hit.iterrows():
                if j <= i:
                    continue
                sim = prefix_similarity(hi['ref_name'], hj['ref_name'])
                internal_rows.append({
                    'hit1': hi['ref_name'],
                    'hit2': hj['ref_name'],
                    'similarity': sim
                })
        internal_df = pd.DataFrame(internal_rows)
#        print("Internal prefix similarity:")
#        print(internal_df)

        # EXTERNAL PREFIX SIMILARITY MATRIX
        external_rows = []
        if not opp.empty:
            for _, hi in top_hit.iterrows():
                for _, o in opp.iterrows():
                    sim = prefix_similarity(hi['ref_name'], o['ref_name'])
                    external_rows.append({
                        'top_hit': hi['ref_name'],
                        'opp_hit': o['ref_name'],
                        'similarity': sim
                    })
        external_df = pd.DataFrame(external_rows)
#        print("External prefix similarity:")
#        print(external_df)

        opp_names = opp['ref_name'].tolist()
        keep_opp = []
        removed_opp = []
        for o_name in opp_names:
            sim = max(prefix_similarity(o_name, t) for t in top_names)
            if sim >= 0.7:
                keep_opp.append(o_name)
            else:
                removed_opp.append(o_name)
        if removed_opp:
            print(f"Removed opp hits below threshold: {removed_opp}")
        opp = opp[opp['ref_name'].isin(keep_opp)]

        internal_prefix = None
        internal_len = 0
        for _, row in internal_df.iterrows():
            if row['similarity'] > internal_len:
                internal_len = row['similarity']
                internal_prefix = lcp(row['hit1'], row['hit2'])
        external_prefix = None
        external_len = 0
        for _, row in external_df.iterrows():
            if row['similarity'] > external_len:
                external_len = row['similarity']
                external_prefix = lcp(row['top_hit'], row['opp_hit'])

        # Condition for known gDNA
        top_hit['percent'] = pd.to_numeric(top_hit['percent'], errors='coerce')
        top_hit['vs_ref'] = pd.to_numeric(top_hit['vs_ref'], errors='coerce')

        # Filter only gDNA rows
        gDNA_rows = top_hit[top_hit['library'] == 'gDNA']

        # Now check only those rows
        if not gDNA_rows.empty and (gDNA_rows['percent'] == 100.0).all() and (gDNA_rows['vs_ref'] == 0).all():
            chosen = gDNA_rows.iloc[0]['ref_name']
        else:
            # choose the longer prefix
            if external_prefix and external_len > internal_len:
                chosen = external_prefix
            else:
                chosen = internal_prefix
        gene_name = chosen.rstrip('_- ') if chosen else None
    else:
        gene_name = top_hit.iloc[0]['ref_name']

    dfs_to_concat = [df for df in [top_hit, opp] if not df.empty]
    final_cluster_df = pd.concat(dfs_to_concat, ignore_index=True) if dfs_to_concat else pd.DataFrame()
    avg_pct, start, end = _cluster_metrics(final_cluster_df)
    
    # collect unique sources from top_hit and any filtered_opp rows
    top_sources = set(top_hit['source'].tolist())
    opp_sources = set(opp['source'].tolist()) if 'opp' in locals() else set()
    all_sources = top_sources.union(opp_sources)

    print ("---0---" * 10)
    print (final_cluster_df)
    print ("-*-*-*-" * 10)
    # take first letter of each, lowercase, sort, and join with '+'
    source_code = "+".join(sorted(s[0].lower() for s in all_sources))

    # Ensure correct types for comparison
    final_cluster_df['percent'] = pd.to_numeric(final_cluster_df['percent'], errors='coerce')
    final_cluster_df['vs_ref'] = pd.to_numeric(final_cluster_df['vs_ref'], errors='coerce')

    gDNA_check = final_cluster_df[
        (final_cluster_df['library'] == 'gDNA') &
        (final_cluster_df['percent'] == 100.0) &
        (final_cluster_df['vs_ref'] == 0)
        ]

    if not gDNA_check.empty:
        return gene_name, avg_pct, start, end, source_code
    else:
        return f"{gene_name}_like", avg_pct, start, end, source_code

def _best_hit2(df: pd.DataFrame, exclude_source=None) -> pd.DataFrame:
    """
    Return all ties for the best hit:
      - If any hits with vs_ref == 0 and percent == 100 exist, return all such hits.
      - Otherwise, return all hits with minimal vs_ref and maximal percent (ties allowed).
    Optionally exclude one or more sources (string, list-like, Series, or Index).
    """
    df_copy = df.copy()
    # handle exclusion of sources
    if exclude_source is not None:
        # scalar string
        if isinstance(exclude_source, str):
            df_copy = df_copy[df_copy['source'] != exclude_source]
        else:
            # assume iterable of sources (list, set, tuple, Series, Index)
            try:
                excl_list = list(exclude_source)
            except Exception:
                excl_list = [exclude_source]
            df_copy = df_copy[~df_copy['source'].isin(excl_list)]

    # perfect matches
    perfect = df_copy[(df_copy['vs_ref'] == 0) & (df_copy['percent'] == 100)]
    if not perfect.empty:
        return perfect.copy()

    # find minimal vs_ref
    min_vs = df_copy['vs_ref'].min()
    candidates = df_copy[df_copy['vs_ref'] == min_vs]
    # within those, maximal percent
    max_pct = candidates['percent'].astype(float).max()
    best = candidates[candidates['percent'].astype(float) == max_pct]
    return best.copy()

def process_hits(top_hit: pd.DataFrame, opp_hit: pd.DataFrame):

    # Join into final_cluster_df AFTER all filtering
    dfs_to_concat = [df for df in [top_hit, opp_hit] if not df.empty]
    final_cluster_df = pd.concat(dfs_to_concat, ignore_index=True) if dfs_to_concat else pd.DataFrame()

    # Ensure numeric columns before filtering
    top_hit['percent'] = pd.to_numeric(top_hit['percent'], errors='coerce')
    top_hit['vs_ref'] = pd.to_numeric(top_hit['vs_ref'], errors='coerce')

    # Early return if top_hit contains gDNA hit with perfect match
    gDNA_rows = top_hit[
        (top_hit['library'] == 'gDNA') &
        (top_hit['percent'] == 100.0) &
        (top_hit['vs_ref'] == 0)
    ]

    avg_pct, start, end = _cluster_metrics(final_cluster_df)
    sources = set(final_cluster_df['source'].tolist())
    source_code = "+".join(sorted(s[0].lower() for s in sources))

    # Condition for perfect gDNA
    if not gDNA_rows.empty:
        gene_name = gDNA_rows.iloc[0]['ref_name']
#        print (gene_name)     
    else:
        ref_names = final_cluster_df['ref_name'].tolist()

        if len(ref_names) == 1:
            gene_name = ref_names[0]
        else:
            gene_name = ref_names[0]
            for other in ref_names[1:]:
                gene_name = lcp(gene_name, other)
            # Clean up 
            gene_name = gene_name.rstrip('_- ') if gene_name else None
        gene_name = gene_name + "_like"
    return gene_name, avg_pct, start, end, source_code


def _best_hit(df: pd.DataFrame, exclude_source=None) -> pd.DataFrame:
    df_copy = df.copy()
    df_copy['percent'] = pd.to_numeric(df_copy['percent'], errors='coerce')
    df_copy['vs_ref'] = pd.to_numeric(df_copy['vs_ref'], errors='coerce')

    # handle exclusion of sources
    if exclude_source is not None:
        # scalar string
        if isinstance(exclude_source, str):
            df_copy = df_copy[df_copy['source'] != exclude_source]
        else:
            # assume iterable of sources (list, set, tuple, Series, Index)
            try:
                excl_list = list(exclude_source)
            except Exception:
                excl_list = [exclude_source]
            df_copy = df_copy[~df_copy['source'].isin(excl_list)]
        return df_copy.copy()

    # top hit detection
    perfects = df_copy[(df_copy['vs_ref'] == 0) & (df_copy['percent'] == 100)]
    if not perfects.empty:
        return perfects.copy()
    else:
        # Pick top percent
        max_percent = df_copy['percent'].max()
        pct_top = df_copy[df_copy['percent'] == max_percent]
        # Step 2: among those, filter by highest align_percent
        max_align = pct_top['align_percent'].max()
        bests = pct_top[pct_top['align_percent'] == max_align]
        return bests.copy()
'''
# Top hit detection
    perfects = df_copy[(df_copy['vs_ref'] == 0) & (df_copy['percent'] == 100)]
    if not perfects.empty:
        bests = perfects.copy()
    else:
        max_percent = df_copy['percent'].max()
        pct_top = df_copy[df_copy['percent'] == max_percent]
        max_align = pct_top['align_percent'].max()
        bests = pct_top[pct_top['align_percent'] == max_align]

    # === Additional logic: if top hit is "protein", prefer alternative if available ===
    if len(bests) == 1 and bests.iloc[0]['source'] == 'protein':
        # Look for any other hits not from protein
        others = df_copy[df_copy['source'] != 'protein']
        if not others.empty:
            # Recursively re-call function excluding 'protein'
            return _best_hit(df, exclude_source='protein')

    return bests.copy()
'''

def get_common_name2(name1: str, name2: str) -> str:
    min_len = min(len(name1), len(name2))
    common = []

    for i in range(min_len):
        if name1[i] == name2[i]:
            common.append(name1[i])
        else:
            break

    # Now ensure we don't cut partial numbers: trim to last digit boundary
    common_str = ''.join(common)

    # Special logic: If there's a '*' in the string, try to include full digits after it
    if '*' in common_str:
        star_idx = common_str.index('*')
        # Keep everything up to star + following digits
        suffix = []
        for c in common_str[star_idx + 1:]:
            if c.isdigit():
                suffix.append(c)
            else:
                break
        return common_str[:star_idx + 1] + ''.join(suffix)
    else:
        return common_str

def get_common_name(name1: str, name2: str) -> str:
    min_len = min(len(name1), len(name2))
    common = []

    for i in range(min_len):
        if name1[i] == name2[i]:
            common.append(name1[i])
        else:
            break

    common_str = ''.join(common)

    # If the names are fully equal
    if name1 == name2:
        return name1

    # If there's a '*', extend past it to include full numeric block
    if '*' in common_str:
        star_idx = common_str.index('*')
        suffix1 = name1[star_idx + 1:]
        suffix2 = name2[star_idx + 1:]
        suffix = []

        for c1, c2 in zip(suffix1, suffix2):
            if c1 == c2 and (c1.isdigit() or c1 == ':'):
                suffix.append(c1)
            else:
                break
        return common_str[:star_idx + 1] + ''.join(suffix)

    return common_str

def count_suffix_after_star(name: str) -> int:
    if '*' in name:
        suffix = name.split('*', 1)[1]
        return len(suffix)
    return 0

def process_cluster(cluster_df: pd.DataFrame):
    # BEST OVERALL HIT
    top_hit = _best_hit(cluster_df)
    best_name, best_src = top_hit['ref_name'], top_hit['source']
    sources = set(cluster_df['source'])

    opp_hit = None

    # MULTI-SOURCE CASE
    opp_hit = _best_hit(cluster_df, exclude_source=best_src)

#    print (top_hit)
#    print (opp_hit)
    # INTERNAL PREFIX SIMILARITY
    internal_rows = []
    top_names = top_hit['ref_name'].tolist()
    for i, hi in top_hit.iterrows():
        for j, hj in top_hit.iterrows():
            if j <= i:
                continue
            common_name = get_common_name(hi['ref_name'], hj['ref_name'])
            internal_rows.append({
                'hit1': hi['ref_name'],
                'hit2': hj['ref_name'],
                'common_name': common_name
            })
    internal_df = pd.DataFrame(internal_rows)
    if not internal_df.empty:
        internal_df['suffix_len'] = internal_df['common_name'].apply(count_suffix_after_star)

    # EXTERNAL PREFIX SIMILARITY
    external_rows = []
    if not opp_hit.empty:
        for _, hi in top_hit.iterrows():
            for _, o in opp_hit.iterrows():
                common_name = get_common_name(hi['ref_name'], o['ref_name'])
                external_rows.append({
                    'top_hit': hi['ref_name'],
                    'opp_hit': o['ref_name'],
                    'common_name': common_name
                })
    external_df = pd.DataFrame(external_rows)
    if not external_df.empty:
        external_df['suffix_len'] = external_df['common_name'].apply(count_suffix_after_star)

    removed_opp = []
    
#    print (internal_df)
#    print (external_df)

    if not opp_hit.empty:
        duplicate_libraries = opp_hit['library'].value_counts()
        dup_libs = duplicate_libraries[duplicate_libraries > 1].index.tolist()

        for lib in dup_libs:
            subset = opp_hit[opp_hit['library'] == lib]
            best_ref = None
            best_sim = -1
            for _, row in subset.iterrows():
                sim = max(prefix_similarity(row['ref_name'], t) for t in top_names)
                if sim > best_sim:
                    best_sim = sim
                    best_ref = row['ref_name']
            # Identify and track rows to remove
            to_remove = subset[subset['ref_name'] != best_ref]
            removed_opp.extend(to_remove['ref_name'].tolist())

            # Remove them from opp_hit
            opp_hit = opp_hit[~((opp_hit['library'] == lib) & (opp_hit['ref_name'] != best_ref))]

        if not external_df.empty:
            low_suffix_refs = external_df[external_df['suffix_len'] < 3]['opp_hit'].unique()
            if len(low_suffix_refs) > 0:
                opp_hit = opp_hit[~opp_hit['ref_name'].isin(low_suffix_refs)]
                removed_opp.extend(low_suffix_refs.tolist())
    
    if removed_opp:
        print(f"Removed opp hits: {removed_opp}")

    return process_hits(top_hit, opp_hit)



# === AGGREGATION ACROSS ALL CLUSTERS ===
def run_all_three(df_cdna: pd.DataFrame, df_gdna: pd.DataFrame, df_protein: pd.DataFrame, out_path: str):
    print("Running analysis on CDNA + GDNA + PROTEIN …")
    cdna_df = df_cdna.copy(); cdna_df['source'] = 'cDNA'
    gdna_df = df_gdna.copy(); gdna_df['source'] = 'gDNA'
    protein_df = df_protein.copy(); protein_df['source'] = 'protein'

    combined = pd.concat([cdna_df, gdna_df, protein_df], ignore_index=True)
    for col in ['contig','strand','roi_start','roi_end','ref_name','percent','source']:
        if col not in combined.columns:
            sys.exit(f"Missing column: {col}")

    clustered = cluster_intervals(combined, 'roi_start', 'roi_end', ['contig','strand'])
    clustered = clustered.sort_values('cluster_id')
    print(clustered.to_string())

    consensus = []

    for cid, dfc in clustered.groupby('cluster_id'):
        strand = dfc.iloc[0]['strand']
        contig = dfc.iloc[0]['contig']
        name, avg_pct, s_coord, e_coord, score = process_cluster(dfc)

        if not name:
            continue  # skip if name is None or empty

        match = re.match(r'^.+\*([a-zA-Z0-9]+)_like', name)
        if match:
            suffix = match.group(1)
            if len(suffix) < 3:
                # special case: rerun _best_hit and override name-based result
                top_hit = _best_hit(dfc)
                top_hit['percent'] = pd.to_numeric(top_hit['percent'], errors='coerce')
                sources = set(top_hit['source'].tolist())
                source_code = "+".join(sorted(s[0].lower() for s in sources))
                for _, row in top_hit.iterrows():
                    gene_name = f"{row['ref_name']}_like_{row['percent']:.2f}"
                    consensus.append({
                        'contig': contig,
                        'roi_start': int(row['roi_start']),
                        'roi_end': int(row['roi_end']),
                        'ref_name': gene_name,
                        'status': source_code,
                        'strand': strand
                    })
                continue  # skip default append, since we already added custom hits

        # fallback: default from process_cluster
        gene_name = f"{name}_{avg_pct:.2f}"
        consensus.append({
            'contig': contig,
            'roi_start': s_coord,
            'roi_end': e_coord,
            'ref_name': gene_name,
            'status': score,
            'strand': strand
        })
    df_cons = pd.DataFrame(consensus)

    if not df_cons.empty:
        df_cons = df_cons.sort_values('roi_start').reset_index(drop=True)

    result = df_cons.sort_values('roi_start')
    print (result)
    result.to_csv(out_path, sep='\t', index=False)

    return result


def _cluster_metrics(h1, h2=None):
    # h1, h2: DataFrame or Series
    if isinstance(h1, pd.Series):
        pcts = [float(h1['percent'])]
        starts = [h1['roi_start']]
        ends = [h1['roi_end']]
    else:
        pcts = h1['percent'].astype(float).tolist()
        starts = h1['roi_start'].tolist()
        ends = h1['roi_end'].tolist()
    if h2 is not None:
        if isinstance(h2, pd.Series):
            pcts += [float(h2['percent'])]
            starts += [h2['roi_start']]
            ends += [h2['roi_end']]
        else:
            pcts += h2['percent'].astype(float).tolist()
            starts += h2['roi_start'].tolist()
            ends += h2['roi_end'].tolist()
    avg_pct = sum(pcts) / len(pcts)
    return avg_pct, int(min(starts)), int(max(ends))


def process_cluster_protein(cluster_df):
    print(cluster_df.to_string())

    # Step 1: find all tied “best” hits, excluding 'protein'
    best_hits = _best_hit(cluster_df, exclude_source='protein')
    # If multiple best hits share the same source, return all names as a list,
    # with the same metrics and source.
    if len(best_hits) > 1:
        sources = set(best_hits['source'])
        if len(sources) == 1:
            # All tied hits are from the same source
            source = next(iter(sources))
            # Convert vs_ref and percent to numeric on the first row (they're identical for all ties)
            first = best_hits.iloc[0]
            first['vs_ref']  = pd.to_numeric(first['vs_ref'],  errors='raise')
            first['percent'] = pd.to_numeric(first['percent'], errors='raise')

            # Compute metrics once
            avg_pct, start, end = _cluster_metrics(first)

            # Collect all tied ref_names
            names = list(best_hits['ref_name'])
            return names, avg_pct, start, end, source

    # Otherwise, pick exactly one “best” hit (the first tied row if there are multiple sources)
    best = best_hits.iloc[0]
    best['vs_ref']  = pd.to_numeric(best['vs_ref'],  errors='raise')
    best['percent'] = pd.to_numeric(best['percent'], errors='raise')

    best_name, best_src = best['ref_name'], best['source']
    sources = set(cluster_df['source'])
    print(best_name, best_src)

    # If there is only one source in the entire cluster, return that single hit’s metrics
    if len(sources) == 1:
        avg_pct, start, end = _cluster_metrics(best)
        return best_name, avg_pct, start, end, best_src

    # Both cDNA and gDNA present → find the best opposite-source hit
    opp_hits = _best_hit(cluster_df, exclude_source=best_src)
    opp = opp_hits.iloc[0]  # pick the first if there is a tie among opposite-source hits

    avg_pct, start, end = _cluster_metrics(best, opp)

    score = "gt"
    if best_src == 'gDNA' and best['vs_ref'] == 0 and best['percent'] == 100.00:
        # identical names → use it directly
        if best_name == opp['ref_name']:
            return best_name, avg_pct, start, end, score

        # different names → commonprefix + '_like'
        prefix = commonprefix([best_name, opp['ref_name']]).rstrip('_- ')
        return f"{prefix}_like", avg_pct, start, end, score

    # All other mixed cases → “best_name_like”
    return f"{best_name}_like", avg_pct, start, end, score
 

def cluster_intervals2(df, start_col, end_col, group_cols):
    """
    Cluster intervals based on gDNA if present; otherwise cluster all rows.
    Removes cDNA/protein hits spanning >1 gDNA cluster.
    """
    df_copy = df.copy()
    df_copy['cluster_id'] = None

    # If gDNA present, cluster gDNA only
    if 'gDNA' in df_copy['source'].unique():
        gdf = df_copy[df_copy['source'] == 'gDNA'].copy()
        gdf['cluster_id'] = None
        # cluster gDNA
        for grp_vals, grp_df in gdf.groupby(group_cols):
            sorted_idx = grp_df[start_col].sort_values().index
            clusters, current, current_end = [], [], None
            for idx in sorted_idx:
                s, e = grp_df.at[idx, start_col], grp_df.at[idx, end_col]
                if current_end is None or s > current_end:
                    if current: clusters.append(current)
                    current, current_end = [idx], e
                else:
                    current.append(idx)
                    current_end = max(current_end, e)
            if current: clusters.append(current)

            for cid, cluster in enumerate(clusters, start=1):
                label = f"cluster_{'_'.join(map(str, grp_vals))}_{cid}"
                for idx in cluster:
                    gdf.at[idx, 'cluster_id'] = label

        # assign back
        df_copy.loc[gdf.index, 'cluster_id'] = gdf['cluster_id']

        # drop multi-cluster spans
        drop_idx = []
        for idx, row in df_copy[df_copy['source'] != 'gDNA'].iterrows():
            mask = True
            for col in group_cols:
                mask &= (gdf[col] == row[col])
            overlaps = gdf[mask & ~((row[end_col] < gdf[start_col]) |
                                     (row[start_col] > gdf[end_col]))]
            cids = overlaps['cluster_id'].dropna().unique()
            if len(cids) > 1:
                print(f"Removing {row['source']} idx={idx}, spans {list(cids)}")
                drop_idx.append(idx)
            elif len(cids) == 1:
                df_copy.at[idx, 'cluster_id'] = cids[0]
        df_copy.drop(index=drop_idx, inplace=True)
        return df_copy

    # no gDNA: cluster all
    # Keep this for now. Have to see if any case where protein and cDNA only which span over each other.
    for grp_vals, grp_df in df_copy.groupby(group_cols):
        sorted_idx = grp_df[start_col].sort_values().index
        clusters, current, current_end = [], [], None
        for idx in sorted_idx:
            s, e = grp_df.at[idx, start_col], grp_df.at[idx, end_col]
            if current_end is None or s > current_end:
                if current: clusters.append(current)
                current, current_end = [idx], e
            else:
                current.append(idx)
                current_end = max(current_end, e)
        if current: clusters.append(current)
        for cid, cluster in enumerate(clusters, start=1):
            label = f"cluster_{'_'.join(map(str, grp_vals))}_{cid}"
            for idx in cluster:
                df_copy.at[idx, 'cluster_id'] = label
    return df_copy


def cluster_intervals3(df, start_col, end_col, group_cols):
    """
    Cluster intervals based on gDNA if present; otherwise cluster all rows.
    Removes cDNA/protein hits spanning >1 gDNA cluster.
    Clusters remaining unassigned intervals as fallback.
    """
    df_copy = df.copy()
    df_copy['cluster_id'] = None

    # If gDNA present, cluster gDNA only
    if 'gDNA' in df_copy['source'].unique():
        gdf = df_copy[df_copy['source'] == 'gDNA'].copy()
        gdf['cluster_id'] = None

        for grp_vals, grp_df in gdf.groupby(group_cols):
            sorted_idx = grp_df[start_col].sort_values().index
            clusters, current, current_end = [], [], None
            for idx in sorted_idx:
                s, e = grp_df.at[idx, start_col], grp_df.at[idx, end_col]
                if current_end is None or s > current_end:
                    if current: clusters.append(current)
                    current, current_end = [idx], e
                else:
                    current.append(idx)
                    current_end = max(current_end, e)
            if current: clusters.append(current)

            for cid, cluster in enumerate(clusters, start=1):
                label = f"cluster_{'_'.join(map(str, grp_vals))}_{cid}"
                for idx in cluster:
                    gdf.at[idx, 'cluster_id'] = label

        df_copy.loc[gdf.index, 'cluster_id'] = gdf['cluster_id']

        # drop multi-cluster spans
        drop_idx = []
        for idx, row in df_copy[df_copy['source'] != 'gDNA'].iterrows():
            mask = True
            for col in group_cols:
                mask &= (gdf[col] == row[col])
            overlaps = gdf[mask & ~((row[end_col] < gdf[start_col]) |
                                     (row[start_col] > gdf[end_col]))]
            cids = overlaps['cluster_id'].dropna().unique()
            if len(cids) > 1:
                print(f"Removing {row['source']} idx={idx}, spans {list(cids)}")
                drop_idx.append(idx)
            elif len(cids) == 1:
                df_copy.at[idx, 'cluster_id'] = cids[0]
        df_copy.drop(index=drop_idx, inplace=True)

    # Fallback: cluster remaining unassigned intervals
    fallback_df = df_copy[df_copy['cluster_id'].isna()]
    for grp_vals, grp_df in fallback_df.groupby(group_cols):
        sorted_idx = grp_df[start_col].sort_values().index
        clusters, current, current_end = [], [], None
        for idx in sorted_idx:
            s, e = grp_df.at[idx, start_col], grp_df.at[idx, end_col]
            if current_end is None or s > current_end:
                if current: clusters.append(current)
                current, current_end = [idx], e
            else:
                current.append(idx)
                current_end = max(current_end, e)
        if current: clusters.append(current)

        for cid, cluster in enumerate(clusters, start=1):
            label = f"fallback_cluster_{'_'.join(map(str, grp_vals))}_{cid}"
            for idx in cluster:
                df_copy.at[idx, 'cluster_id'] = label

    return df_copy

def cluster_intervals4(df, start_col, end_col, group_cols):
    """
    Cluster intervals based on gDNA if present; otherwise cluster all rows.
    Removes cDNA/protein hits spanning >1 gDNA cluster.
    Clusters remaining unassigned intervals using overlap-aware fallback.
    """
    df_copy = df.copy()
    df_copy['cluster_id'] = None

    # Ensure start/end are integers
    df_copy[start_col] = df_copy[start_col].astype(int)
    df_copy[end_col] = df_copy[end_col].astype(int)

    # Helper function for interval clustering
    def build_clusters(group_df, label_prefix):
        sorted_df = group_df.sort_values(by=start_col)
        clusters = []
        current, current_end = [], None
        for idx, row in sorted_df.iterrows():
            s, e = row[start_col], row[end_col]
            if current_end is None or s > current_end:
                if current:
                    clusters.append(current)
                current = [idx]
                current_end = e
            else:
                current.append(idx)
                current_end = max(current_end, e)
        if current:
            clusters.append(current)
        # Assign cluster ids
        cluster_labels = {}
        for cid, cluster in enumerate(clusters, start=1):
            label = f"{label_prefix}_{cid}"
            for idx in cluster:
                cluster_labels[idx] = label
        return cluster_labels

    # Step 1: Cluster gDNA (if present)
    if 'gDNA' in df_copy['source'].unique():
        gdf = df_copy[df_copy['source'] == 'gDNA'].copy()
        g_labels = {}
        for grp_vals, grp_df in gdf.groupby(group_cols):
            prefix = f"cluster_{'_'.join(map(str, grp_vals))}"
            g_labels.update(build_clusters(grp_df, prefix))
        df_copy.loc[g_labels.keys(), 'cluster_id'] = pd.Series(g_labels)

        # Step 2: Assign non-gDNA rows based on overlap
        non_gdf = df_copy[df_copy['source'] != 'gDNA']
        gdf = df_copy[df_copy['source'] == 'gDNA'].copy()  # clustered

        drop_idx = []
        for idx, row in non_gdf.iterrows():
            mask = True
            for col in group_cols:
                mask &= (gdf[col] == row[col])
            overlaps = gdf[mask & ~((row[end_col] < gdf[start_col]) |
                                     (row[start_col] > gdf[end_col]))]
            cids = overlaps['cluster_id'].dropna().unique()
            if len(cids) > 1:
                print(f"Removing {row['source']} idx={idx}, spans {list(cids)}")
                drop_idx.append(idx)
            elif len(cids) == 1:
                df_copy.at[idx, 'cluster_id'] = cids[0]

        df_copy.drop(index=drop_idx, inplace=True)

    # Step 3: Fallback clustering for still-unassigned rows
    fallback_df = df_copy[df_copy['cluster_id'].isna()]
    f_labels = {}
    for grp_vals, grp_df in fallback_df.groupby(group_cols):
        prefix = f"fallback_cluster_{'_'.join(map(str, grp_vals))}"
        f_labels.update(build_clusters(grp_df, prefix))
    df_copy.loc[f_labels.keys(), 'cluster_id'] = pd.Series(f_labels)

    return df_copy

def cluster_intervals(df, start_col, end_col, group_cols):
    """
    Cluster gDNA if available, assign overlapping cDNA/protein hits,
    and chain unclustered hits that overlap clustered transcripts.
    """
    import pandas as pd

    df_copy = df.copy()
    df_copy['cluster_id'] = None
    df_copy[start_col] = df_copy[start_col].astype(int)
    df_copy[end_col] = df_copy[end_col].astype(int)

    def build_clusters(group_df, label_prefix):
        sorted_df = group_df.sort_values(by=start_col)
        clusters = []
        current, current_end = [], None
        for idx, row in sorted_df.iterrows():
            s, e = row[start_col], row[end_col]
            if current_end is None or s > current_end:
                if current:
                    clusters.append(current)
                current = [idx]
                current_end = e
            else:
                current.append(idx)
                current_end = max(current_end, e)
        if current:
            clusters.append(current)
        cluster_labels = {}
        for cid, cluster in enumerate(clusters, start=1):
            label = f"{label_prefix}_{cid}"
            for idx in cluster:
                cluster_labels[idx] = label
        return cluster_labels

    # Step 1: Cluster gDNA
    if 'gDNA' in df_copy['source'].unique():
        gdf = df_copy[df_copy['source'] == 'gDNA']
        g_labels = {}
        for grp_vals, grp_df in gdf.groupby(group_cols):
            prefix = f"cluster_{'_'.join(map(str, grp_vals))}"
            g_labels.update(build_clusters(grp_df, prefix))
        df_copy.loc[g_labels.keys(), 'cluster_id'] = pd.Series(g_labels)

        # Step 2: Assign cDNA/protein by overlap to gDNA
        non_gdf = df_copy[df_copy['source'] != 'gDNA']
        gdf = df_copy[df_copy['source'] == 'gDNA']

        drop_idx = []
        for idx, row in non_gdf.iterrows():
            mask = True
            for col in group_cols:
                mask &= (gdf[col] == row[col])
            overlaps = gdf[mask & ~((row[end_col] < gdf[start_col]) |
                                     (row[start_col] > gdf[end_col]))]
            cids = overlaps['cluster_id'].dropna().unique()
            if len(cids) > 1:
                drop_idx.append(idx)
            elif len(cids) == 1:
                df_copy.at[idx, 'cluster_id'] = cids[0]
        df_copy.drop(index=drop_idx, inplace=True)

    # Step 3: Assign unclustered hits based on overlap with already-clustered rows (chaining)
    unclustered = df_copy[df_copy['cluster_id'].isna()]
    clustered = df_copy[df_copy['cluster_id'].notna()]

    for idx, row in unclustered.iterrows():
        mask = True
        for col in group_cols:
            mask &= (clustered[col] == row[col])
        overlaps = clustered[mask & ~((row[end_col] < clustered[start_col]) |
                                      (row[start_col] > clustered[end_col]))]
        cids = overlaps['cluster_id'].dropna().unique()
        if len(cids) == 1:
            df_copy.at[idx, 'cluster_id'] = cids[0]

    # Step 4: Fallback cluster remaining unclustered entries
    final_unclustered = df_copy[df_copy['cluster_id'].isna()]
    f_labels = {}
    for grp_vals, grp_df in final_unclustered.groupby(group_cols):
        prefix = f"fallback_cluster_{'_'.join(map(str, grp_vals))}"
        f_labels.update(build_clusters(grp_df, prefix))
    df_copy.loc[f_labels.keys(), 'cluster_id'] = pd.Series(f_labels)

    return df_copy


def run_cdna_gdna(df_cdna, df_gdna, out_path):
    print("Running analysis on CDNA + GDNA …")

    cdna_df = df_cdna.copy(); cdna_df['source'] = 'cDNA'
    gdna_df = df_gdna.copy(); gdna_df['source'] = 'gDNA'
    combined = pd.concat([cdna_df, gdna_df], ignore_index=True)
#    print (combined)
    for col in ['contig','strand','roi_start','roi_end','ref_name','percent','source']:
        if col not in combined.columns: sys.exit(f"Missing column: {col}")
    clustered = cluster_intervals(combined,'roi_start','roi_end',['contig','strand'])

#    print (clustered)
    consensus = []
    for cid, dfc in clustered.groupby('cluster_id'):
        # extract contig and strand
        strand = dfc.iloc[0]['strand']
        contig = dfc.iloc[0]['contig']

        name, avg_pct, s_coord, e_coord, score = process_cluster(dfc)
        if name:
            gene_name = f"{name}_{avg_pct:.2f}"
            consensus.append({
                'contig': contig,
                'gene_name': gene_name,
                'start': s_coord,
                'end': e_coord,
                'score': score,
                'strand': strand,
            })
    # create DataFrame with merged gene_name; drop original name & avg_percent
    df_cons = pd.DataFrame(consensus)
    df_cons = df_cons.sort_values("start")
    print("\nConsensus summary:")
    print(df_cons.to_string(index=False))

    if out_path:
        df_cons.to_csv(f"{out_path}", sep='\t', index=False, header=False)
        print(f"Saved consensus to {out_path}")


def lcp2(a: str, b: str) -> str:
    """Longest common prefix."""
    i = 0
    for ca, cb in zip(a, b):
        if ca == cb:
            i += 1
        else:
            break
    return a[:i]


def base_name(s: str) -> str:
    """Prefix before the first '*'"""
    return s.split('*', 1)[0]


def _select_best(group: pd.DataFrame, src: str) -> pd.DataFrame:
    """Select best hits for a given source by vs_ref then percent."""
    rows = group[group['source'] == src]
    min_vs = rows['vs_ref'].min()
    candidates = rows[rows['vs_ref'] == min_vs]
    max_pct = candidates['percent'].astype(float).max()
    return candidates[candidates['percent'].astype(float) == max_pct]


def _merge_hits(dna_row, prot_row, label: str) -> dict:
    """Merge two rows (DNA & protein) into a final annotation."""
    # choose winner for naming
    winner = dna_row if float(dna_row['percent']) >= float(prot_row['percent']) else prot_row
    pct = float(winner['percent'])
    name = winner['ref_name']
    suffix = f"_{pct:.0f}p"
    gene_name = f"{name}{'_like' if 'novel' in (dna_row['status'], prot_row['status']) else ''}{suffix}"

    contig = winner['contig']
    start = min(dna_row['roi_start'], prot_row['roi_start'])
    end = max(dna_row['roi_end'], prot_row['roi_end'])
    strand = prot_row['strand'] if prot_row['strand'] != dna_row['strand'] else dna_row['strand']
    score = f"{label}+p"

    return dict(
        contig=contig,
        roi_start=start,
        roi_end=end,
        gene_name=gene_name,
        score=score,
        strand=strand
    )

def _single_source_hit(row, label: str) -> dict:
    """Create gene record for single-source clusters."""
    base = row['ref_name']
    pct = float(row['percent'])
    if row['status'] == 'known':
        suffix = '_100' if label == 'p' else '_100p'
    else:
        suffix = f"_like_{pct:.0f}p"

    return dict(
        contig=row['contig'],
        roi_start=row['roi_start'],
        roi_end=row['roi_end'],
        gene_name=f"{base}{suffix}",
        score=label,
        strand=row['strand']
    )

def run_dna_protein(df_dna: pd.DataFrame,
                    df_protein: pd.DataFrame,
                    out_path: str,
                    dna_label: str = 'cDNA') -> pd.DataFrame:
    """
    Generic handler for DNA (cDNA/gDNA) vs protein analysis.
    dna_label: either 'cDNA' or 'gDNA'
    writes a TSV to out_path and returns the DataFrame.
    """
    print(f"Running analysis on {dna_label} + PROTEIN …")
    df_dna = df_dna.copy()
    df_protein = df_protein.copy()
    df_dna['source'] = dna_label
    df_protein['source'] = 'protein'

    # 1) Cluster DNA intervals
    clusters = cluster_intervals(df_dna, 'roi_start', 'roi_end', ['contig', 'strand'])
    cid = f"{dna_label[0].lower()}_cluster_id"
    clusters = clusters.rename(columns={'cluster_id': cid})

    # 2) Assign protein to clusters
    prot = df_protein.copy()
    prot[cid] = None
    for idx, prow in prot.iterrows():
        same = clusters[(clusters['contig'] == prow['contig']) &
                        (clusters['strand'] == prow['strand'])]
        overlaps = same[~((prow['roi_end'] < same['roi_start']) |
                          (prow['roi_start'] > same['roi_end']))]
        if len(overlaps) == 1:
            prot.at[idx, cid] = overlaps.iloc[0][cid]
    prot = prot.dropna(subset=[cid])

    # 3) Combine and sort
    combined = pd.concat([clusters, prot], ignore_index=True, sort=False)
    combined = combined.sort_values(by=cid)

    # 4) Best-hit logic
    records = []
    for cluster_id, group in combined.groupby(cid):
        sources = set(group['source'])
        label = dna_label[0].lower()  # 'c' or 'g'

        if {'protein', dna_label}.issubset(sources):
            # both sources
            best_dna = _select_best(group, dna_label)
            best_prot = _select_best(group, 'protein')

            # try find matching base_name
            pair = None
            for _, d in best_dna.iterrows():
                for _, p in best_prot.iterrows():
                    if base_name(d['ref_name']) == base_name(p['ref_name']):
                        pair = (d, p)
                        break
                if pair:
                    break

            if pair:
                rec = _merge_hits(pair[0], pair[1], label)
            else:
                # fallback to DNA-only
                rec = _single_source_hit(best_dna.iloc[0], label)
            records.append(rec)
            continue

        if sources == {'protein'}:
            row = group[group['source'] == 'protein']
            row = row.loc[row['ref_name'].str.len().idxmax()]
            records.append(_single_source_hit(row, 'p'))
            continue

        if sources == {dna_label}:
            row = group[group['source'] == dna_label]
            row = row.loc[row['ref_name'].str.len().idxmax()]
            records.append(_single_source_hit(row, label))
            continue

    result = pd.DataFrame(records,
                           columns=['contig', 'roi_start', 'roi_end', 'gene_name', 'score', 'strand'])
    result = result.sort_values('roi_start')
    result.to_csv(out_path, sep='\t', index=False)
    print(f"  → wrote merged {dna_label}+PROTEIN to '{out_path}'")

    return result


def run_cdna_protein(df_cdna, df_protein, out_path):
    return run_dna_protein(df_cdna, df_protein, out_path, dna_label='cDNA')

def run_gdna_protein(df_gdna, df_protein, out_path):
    return run_dna_protein(df_gdna, df_protein, out_path, dna_label='gDNA')

def run_single_cdna(df_cdna, out_path):
    print("Running analysis on only CDNA …")
    df_cdna["score"] = "cDNA"
    df_cdna=df_cdna[["contig","roi_start","roi_end","gene_name","score","strand"]]
    df_cdna.to_csv(out_path, sep='\t', header=False, index=False)
    print(f"  → output written to {out_path!r}\n")

def run_single_gdna(df_gdna, out_path):
    print("Running analysis on only GDNA …")
    df_gdna["score"] = "gDNA"
    df_gdna=df_gdna[["contig","roi_start","roi_end","gene_name","score","strand"]]
    df_gdna.to_csv(out_path, sep='\t', header=False, index=False)
    print(f"  → output written to {out_path!r}\n")

def run_single_protein(df_protein, out_path):
    print("Running analysis on only PROTEIN …")
    df_protein["score"] = "protein"
    df_protein=df_protein[["contig","roi_start","roi_end","gene_name","score","strand"]]
    df_protein.to_csv(out_path, sep='\t', header=False, index=False)
    print(f"  → output written to {out_path!r}\n")

if __name__ == '__main__':
    main()



