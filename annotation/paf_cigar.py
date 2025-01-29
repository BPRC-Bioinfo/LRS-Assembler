# v0.3
# Giang Le & Jaimy

import pandas as pd
import re
import argparse

# Function to extract matching regions from the CIGAR string and count matches, insertions, and deletions
def extract_match_regions_and_counts(cigar: str, query_start: int, target_start: int):
    cigar_elements = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    query_pos = query_start
    target_pos = target_start
    # Lists to store the extracted data
    query_starts = []
    query_ends = []
    target_starts = []
    target_ends = []
    match_counts = []
    insertion_counts = []
    deletion_counts = []

    # Track merge condition
    merged_query_start = None
    merged_target_start = None
    is_merging = False  # Flag to track merging state
    match_count = 0
    insertion_count = 0
    deletion_count = 0

    for length, operation in cigar_elements:
        length = int(length)

        # If the operation is M/I/D and we are not merging, start merging
        if operation in 'MID':
            if not is_merging:
                merged_query_start = query_pos
                merged_target_start = target_pos
                is_merging = True

            # Update query and target positions based on the operation
            if operation == 'M':
                query_pos += length
                target_pos += length
                match_count += length
            elif operation == 'I':
                query_pos += length
                insertion_count += length
            elif operation == 'D':
                target_pos += length
                deletion_count += length

        # If the operation is N (a gap/skip), we finalize the merge and start a new region
        elif operation == 'N':
            if is_merging:
                # Store the merged region
                query_starts.append(merged_query_start)
                query_ends.append(query_pos)
                target_starts.append(merged_target_start)
                target_ends.append(target_pos)
                match_counts.append(match_count)
                insertion_counts.append(insertion_count)
                deletion_counts.append(deletion_count)

                # Reset the merging state and counters
                is_merging = False
                match_count = 0
                insertion_count = 0
                deletion_count = 0

            # Skip N operations and don't store them
            target_pos += length

    # Finalize any remaining merge
    if is_merging:
        query_starts.append(merged_query_start)
        query_ends.append(query_pos)
        target_starts.append(merged_target_start)
        target_ends.append(target_pos)
        match_counts.append(match_count)
        insertion_counts.append(insertion_count)
        deletion_counts.append(deletion_count)

    # Return a list of rows instead of a DataFrame
    return [
        {
            'Query_start': qs,
            'Query_end': qe,
            'Target_start': ts,
            'Target_end': te,
            'Matches': mc,
            'Insertions': ic,
            'Deletions': dc
        }
        for qs, qe, ts, te, mc, ic, dc in zip(query_starts, query_ends, target_starts, target_ends, match_counts, insertion_counts, deletion_counts)
    ]

# Function to process the PAF file and extract match regions along with counts of matches, insertions, and deletions
def process_paf_file(paf_file_path: str) -> list:
    all_rows = []

    with open(paf_file_path, 'r') as paf_file:
        for line in paf_file:
            fields = line.strip().split('\t')

            qname = fields[0]  # Query name (read ID)
            query_start = int(fields[2])  # Query start position
            rname = fields[5]  # Reference name (e.g., chromosome or contig)
            ref_start = int(fields[7])  # Reference start position
            cigar = None

            # Find the 'cg:Z:' tag for the CIGAR string
            for field in fields:
                if field.startswith('cg:Z:'):
                    cigar = field[5:]
                    break

            if cigar is None:
                continue  # Skip lines without a CIGAR string

            # Extract match regions and counts
            merged_coords = extract_match_regions_and_counts(cigar, query_start, ref_start)

            # Add 'CIGAR', 'Query Name', and 'Reference Name' to each row
            for row in merged_coords:
                row['CIGAR'] = cigar
                row['Query_name'] = qname
                row['Target_name'] = rname
                all_rows.append(row)

    return all_rows

def main(input_file: str = None, output_file: str = None):
    parser = argparse.ArgumentParser(description="Extract match regions and count M/I/D from PAF file")
    parser.add_argument("input_file", help="Path to the PAF file")
    parser.add_argument("output_file", help="Path to the output CSV file")

    args = parser.parse_args()

    all_rows = process_paf_file(args.input_file)
    df = pd.DataFrame(all_rows)

    df = df[['Query_name', 'Target_name', 'Query_start', 'Query_end', 'Target_start', 'Target_end', 'CIGAR', 'Matches', 'Insertions', 'Deletions']]
    df.to_csv(args.output_file, index=False, sep = '\t')

if __name__ == "__main__":
    main()
