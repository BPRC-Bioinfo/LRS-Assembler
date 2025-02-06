import pandas as pd
import argparse
import numpy as np

pd.set_option('display.max_columns', None) 


def process_data(info_df, blast_df):
  """
  Processes the given info and blast DataFrames to extract and return relevant data.

  Args:
    info_df: pandas DataFrame containing info data.
    blast_df: pandas DataFrame containing blast data.

  Returns:
    pandas DataFrame containing the processed results.
  """
#  print (blast_df)
  def process_row(row):
    ref, chr, ref_start, ref_end, chr_start, chr_end, cigar, matches, insertions, deletions, strand = row.values
    chr_coords = f"{chr}:{chr_start}-{chr_end}"
    filtered_df = blast_df[(blast_df[0] == chr_coords) & (blast_df[1] == ref) ] 
#    print (filtered_df)
#    print (strand)

    hits = len(filtered_df)
#    print (hits)

    if hits > 0:
      vals = filtered_df.iloc[0, 2:].tolist()
#      print (vals)
    else:
      vals = []

    return pd.Series([ref, ref_start, ref_end, chr_start, chr_end, cigar, matches, insertions, deletions, strand, hits] + vals)

  result_df = info_df.apply(process_row, axis=1)
  
  # Create a list of column names (assuming they are strings)
#  column_names = [i for i in range(12, 23)]
#  print (result_df)
  column_names = [i for i in range(12, 23) if i not in (19, 20)]  # Columns 12 to 23
  result_df[column_names] = result_df[column_names].fillna(-1).astype(int, errors = "ignore")
  result_df.replace(to_replace=-1, value = None, inplace = True)
#  print (result_df)
  
  return result_df

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Process info and blast files.")
  parser.add_argument("info_file", help="Path to the info file.")
  parser.add_argument("blast_file", help="Path to the blast file.")
  parser.add_argument("-o", "--output", default="output.txt", help="Path to the output file.")
  args = parser.parse_args()

  try:
    info_df = pd.read_csv(args.info_file, sep="\t")
    blast_df = pd.read_csv(args.blast_file, sep="\t", header=None)
  except FileNotFoundError:
    print(f"Error: File not found. Please check the file paths: {args.info_file}, {args.blast_file}")
    exit(1)

  result_df = process_data(info_df, blast_df)
  result_df.to_csv(args.output, sep="\t", index=False, header=False)