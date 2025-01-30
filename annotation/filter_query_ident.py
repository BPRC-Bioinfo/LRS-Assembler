import pandas as pd

def filter_longest_aligned_length(df):
    """
    Filters a DataFrame based on aligned length, identity, and other criteria.

    Args:
        df (pd.DataFrame): Input DataFrame with columns like 'Group_name', 'Aligned Length', 
                           'identity_percent', 'mismatch_sum', 'gap_sum', 'chr_start', 'chr_end'.

    Returns:
        pd.DataFrame: Filtered DataFrame with the best hits based on the specified criteria.
    """
    # Sort the dataframe by identity_percent (descending), then by mismatch_sum and gap_sum (ascending)

    pd.set_option('display.max_rows', None)  # Show all rows
    pd.set_option('display.max_columns', None)  # Show all columns (if necessary)
    pd.set_option('display.width', None)  # Auto-adjust the width to display all columns
    pd.set_option('display.max_colwidth', None) 

    print ("Step1: filter by identity_percent, mismatch_sum and gap_sum")  
    df_sorted = df.sort_values(by=['identity_percent', 'mismatch_sum', 'gap_sum'], ascending=[False, True, True])
    print ("Step1: Group by Group_name") 
    df_best_hits = df_sorted.groupby('Group_name', as_index=False).first()

#    print (df_best_hits)

    df_best_hits = df_best_hits.sort_values(by='chr_start', ascending=True)
    df_best_hits_100 = df_best_hits[df_best_hits["identity_percent"] == 100]
    df_best_hits = df_best_hits[df_best_hits["identity_percent"] != 100]
    print ("Step2")
    print (df_best_hits)

    df_best_hits = (df_best_hits.sort_values(["identity_percent", "mismatch_sum"], ascending=[False, True]).groupby(["chr_start", "chr_end"], as_index=False).first()) 
#    print (df_best_hits)
    df_best_hits = (df_best_hits.sort_values(["chr_start", "aligned_percent"], ascending=[False, False]).groupby("chr_start", as_index=False).first()) 
#    print (df_best_hits)
    df_best_hits = (df_best_hits.sort_values(["chr_end", "aligned_percent"], ascending=[False, False]).groupby("chr_end", as_index=False).first())

#    print (df_best_hits)



    df_final = pd.concat([df_best_hits_100, df_best_hits])

    df_final = df_final[df_final['gap_sum'] == 0]

    df_final = (df_final.sort_values(["identity_percent", "mismatch_sum", "insertions", "deletion"], ascending=[False, True, True, True]).groupby(["chr_start", "chr_end"], as_index=False).first()) 
   
    df_final = df_final.sort_values(by=['chr_start', 'identity_percent'], ascending=[True, False])

    def filter_covered_rows(df):
        filtered_rows = []
        current_max_end = -1
        for idx, row in df.iterrows():
            if row['chr_start'] <= current_max_end:
                prev_row = filtered_rows[-1]
                # Vergelijk niet alleen de percentages, maar ook chr_end als percentages gelijk zijn
                if (row['identity_percent'] > prev_row['identity_percent'] or
                    (row['identity_percent'] == prev_row['identity_percent'] and
                     row['aligned_percent'] == prev_row['aligned_percent'] and
                     row['chr_end'] > prev_row['chr_end'])):
                    # Vervang de vorige rij als de nieuwe betere dekking heeft
                    filtered_rows[-1] = row
                    # Update current_max_end naar de grootste chr_end
                    current_max_end = max(current_max_end, row['chr_end'])
            else:
                # Voeg de rij toe als er geen overlap is
                filtered_rows.append(row)
                current_max_end = row['chr_end']  # Update de grens
        return pd.DataFrame(filtered_rows)

    
    df_final = filter_covered_rows(df_final)
    print("-"*100)

    return df_final.sort_values(["chr_start"])


if __name__ == "__main__":
    import argparse

    # Set up argparse to handle command-line arguments
    parser = argparse.ArgumentParser(description="Filter DataFrame based on aligned length and other criteria.")
    parser.add_argument("input_file", help="Path to the input CSV file")
    parser.add_argument("output_file", help="Path to the output Excel file")
    args = parser.parse_args()

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(args.input_file, sep = '\t') 

    # Call the filtering function
    filtered_df = filter_longest_aligned_length(df)

    print (filtered_df)

    # Write the filtered DataFrame to Excel
    filtered_df.to_excel(args.output_file + ".xlsx", index=False, sheet_name="filtered")
    filtered_df.to_csv(args.output_file + ".csv", index=False)
