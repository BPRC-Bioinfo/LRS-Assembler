# ver 0.0.1
# By: Giang Le & Jaimy
# Process paf with blastn information 

import argparse
import pandas as pd

def merging_correcting_coords(bed_file, filter_file, output_file):
    bed_file = bed_file.sort_values(1, ascending = True)
    bed_file.columns = ["contig", "chr_start", "chr_end", "ref_name", "mapq", "strand"]
    bed_file['ref_name'] = bed_file['ref_name'].str.split("_").str[-1]
    bed_file = bed_file.drop(columns = ['contig','mapq'])
          
    bed_file['blast_percent'] = 'Flanking gene'
    bed_file['blast_align'], bed_file['blast_mismatch'], bed_file['blast_gap'] = 0, 0, 0
    bed_file['ref_len'] = bed_file['chr_end'] - bed_file['chr_start']

    bed_file = bed_file[['ref_name', 'chr_start', 'chr_end', 'blast_percent',  'blast_align',  'blast_mismatch',  'blast_gap', 'strand',  'ref_len']]
    bed_file = bed_file.astype({'chr_start': int, 'chr_end': int})
    roi_start = int(bed_file['chr_start'][0]) -1

    filter_file = filter_file.astype({'roi_start': int, 'roi_end': int})
    filter_file['chr_start'] = filter_file['roi_start'] + roi_start
    filter_file['chr_end'] = filter_file['roi_end'] + roi_start
    
    final_dataframe = pd.concat([bed_file, filter_file])
    final_dataframe = final_dataframe.sort_values('chr_start')

    final_dataframe.to_csv(f"{output_file}", sep='\t', index=False) 
    print (final_dataframe)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Group data by the first column and split into sub-groups based on fragment length and coordinates.")
    parser.add_argument("-b","--bed_file", help="Input file containing the data to process")
    parser.add_argument("-f","--filter_file", help="Input file containing the data to process")
    parser.add_argument("-o","--output", help="Output file name (without extension)")

    args = parser.parse_args()

    try:
        bed_file = pd.read_csv(args.bed_file, header = None, sep = "\t")
        filter_file = pd.read_csv(args.filter_file, sep = "\t")
    except FileNotFoundError:
        print("Error: One or all of the files specified were not found (even after initial existence check).")
        exit(1)
    except pd.errors.ParserError:
        print("Error: Could not parse one or both of the input files. Ensure they are valid TSV files.")
        exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading the input files: {e}")
        exit(1)    

    merging_correcting_coords(bed_file, filter_file, args.output)

