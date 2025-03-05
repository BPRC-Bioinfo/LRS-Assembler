# Filter based on blast results
# v0.0.6
# By Giang Le & Jaimy

import pandas as pd

def filter_longest_aligned_length(df):
    pd.set_option('display.max_rows', None)  # Show all rows
    pd.set_option('display.max_columns', None)  # Show all columns (if necessary)
    pd.set_option('display.width', None)  # Auto-adjust the width to display all columns
    pd.set_option('display.max_colwidth', None)

    # Percentage filter
    df_sorted = df.sort_values(by=['blast_percent', 'blast_mismatch', 'blast_gap'], ascending=[False, True, True])
    df_best_hits = df_sorted.groupby('gene_group', as_index=False).first()
    df_best_hits = df_best_hits.sort_values(by='roi_start', ascending=True)
    
    # Novel and Known split
    df_best_hits_100 = df_best_hits[df_best_hits["blast_percent"] == 100]
    df_best_hits = df_best_hits[df_best_hits["blast_percent"] != 100]
    
    # Blast percent & mismatch
    df_best_hits = (df_best_hits.sort_values(["blast_percent", "blast_mismatch"], ascending=[False, True]).groupby(["roi_start", "roi_end"], as_index=False).first()) 
    # Retain the longest blast_align
    df_best_hits = (df_best_hits.sort_values(["roi_start", "blast_align"], ascending=[False, False]).groupby("roi_start", as_index=False).first()) 
    df_best_hits = (df_best_hits.sort_values(["roi_end", "blast_align"], ascending=[False, False]).groupby("roi_end", as_index=False).first())

    # Filtering with 100
    df_top_hits = pd.concat([df_best_hits_100, df_best_hits])
    df_top_hits = (df_top_hits.sort_values(["blast_percent", "blast_mismatch", "minimap_inserts", "minimap_deletions"], ascending=[False, True, True, True]).groupby(["roi_start", "roi_end"], as_index=False).first()) 
    df_top_hits = df_top_hits.sort_values(by=['roi_start', 'blast_percent'], ascending=[True, False])

    def filter_covered_rows(df):
        filtered_rows = []
        current_max_end = -1
        for idx, row in df.iterrows():
            #print(row['gene_group'], row['roi_start'], current_max_end)

            if row['roi_start'] <= current_max_end:
                prev_row = filtered_rows[-1]
                if row['blast_percent'] >= prev_row['blast_percent'] and row['blast_align'] >= prev_row['blast_align']:

                    if row['roi_start'] > prev_row['roi_end']:
                        filtered_rows.append(row)
                    else:
                        filtered_rows[-1] = row
                    current_max_end = row['roi_end']
                else:
                    continue
            else:
                filtered_rows.append(row)
                current_max_end = row['roi_end']
        return pd.DataFrame(filtered_rows)

    df_final = filter_covered_rows(df_top_hits)

    for _, row in df_final.iterrows():
        found = df_sorted[(df_sorted["gene_group"] != row["gene_group"]) & (df_sorted["roi_start"] == row["roi_start"]) & (df_sorted["roi_end"] == row["roi_end"]) & (df_sorted["blast_percent"] == row["blast_percent"]) & (df_sorted["blast_align"] == row["blast_align"])]
        if not found.empty:
            df_final = pd.concat([df_final, found]).reset_index(drop=True)

    return df_final

def add_suffix_to_short_name(df_final):
    group = df_final.copy()
    if len(group) > 1:
        group = group.sort_values(by=['blast_percent', 'blast_align'], ascending=[False, False]).reset_index(drop=True)
        group['ref_name'] = group['ref_name'] + "_like" + (group.index + 1).astype(str)
    return group


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Filter DataFrame based on aligned length and other criteria.")
    parser.add_argument("-i", "--input", help="Path to the input CSV file")
    parser.add_argument("-o", "--output", help="Path to the output Excel file")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep = '\t') 

    filtered_df = filter_longest_aligned_length(df)
    df_final = add_suffix_to_short_name(filtered_df)
    df_final['ref_name'] = df_final['ref_name'].str.rsplit('|', n=1).str[0]
    df_final = df_final.groupby('ref_name', group_keys=False).apply(add_suffix_to_short_name)
    df_final.sort_values(by=['roi_start'], inplace=True)
    df_final.drop('gene_group', axis=1, inplace=True)
    df_final = df_final[['ref_name', 'roi_start', 'roi_end', 'blast_percent', 'blast_align', 'blast_mismatch', 'blast_gap', 'strand', 'contig', 'ref_len']]
    print("Final genes found!")
    print (df_final)


    df_final.to_csv(args.output, sep='\t' ,index=False)
