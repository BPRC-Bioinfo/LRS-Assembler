#!/usr/bin/env python3

# Script to parse gmap outout
# v0.0.1
# By Giang Le

import sys
import re
import pandas as pd

def process_path_line(line):
    # Regex explanation:
    #   - Path (\d+):                captures the gene/path number.
    #   - query (\d+)\.\.(\d+)         captures query start and end.
    #   - \((\d+)\s+bp\)              captures the full reference length.
    #   - =>\s+genome\s+(\S+):        captures the contig name.
    #   - ([\d,]+)\.\.([\d,]+)         captures the ROI coordinates (numbers may include commas).
    #   - \(([-]?\d+)\s+bp\)          captures the ROI length (with an optional negative sign).
    pattern = (
        r"Path (\d+):\s*query\s+(\d+)\.\.(\d+)\s+\((\d+)\s+bp\)\s+=>\s+genome\s+"
        r"(\S+):([\d,]+)\.\.([\d,]+)\s+\(([-]?\d+)\s+bp\)"
    )
    m = re.search(pattern, line)
    if not m:
        return None

    gene_number = int(m.group(1))
    query_start = int(m.group(2))
    query_end = int(m.group(3))
    contig = m.group(5)
    roi_start = int(m.group(6).replace(",", ""))
    roi_end = int(m.group(7).replace(",", ""))
    
    roi_orientation = m.group(8)
    if roi_orientation.startswith('-'):
        strand = '-'
        roi_len = int(roi_orientation[1:])
        # For negative strand, flip the ROI coordinates.
        roi_start, roi_end = roi_end, roi_start
    else:
        strand = '+'
        roi_len = int(roi_orientation)

    return {
        'contig': contig,
        'ref_start': query_start,
        'ref_end': query_end,
        'roi_start': roi_start,
        'roi_end': roi_end,
        'roi_len': roi_len,
        'strand': strand,
        'gene_number': gene_number
    }

def parse_data(data):

    result = {}
    current_record = None  
    current_path = None 

    for line in data.splitlines():
        if not line.strip():
            continue 

        if line.startswith('>'):
            record_name = line[1:].strip()
            result[record_name] = {}
            current_record = result[record_name]
            current_path = None
            continue

        stripped_line = line.strip()

        if stripped_line.startswith("Path "):
            path_details = process_path_line(stripped_line)
            if path_details and current_record is not None:
                current_record.update(path_details)
            continue

        if (line.startswith("    ") or line.startswith("\t")) and current_record is not None:
            prop_line = line.strip()
            if prop_line.startswith("Number of exons:"):
                try:
                    exon_val = int(prop_line.split("Number of exons:")[1].strip())
                    current_record["exons"] = exon_val
                except ValueError:
                    current_record["exons"] = None
                continue

            # Parse the "Percent identity:" line.
            if prop_line.startswith("Percent identity:"):
                percent_pattern = (
                    r"Percent identity:\s*([\d.]+)\s+\((\d+)\s+matches,\s+"
                    r"(\d+)\s+mismatches,\s+(\d+)\s+indels,\s+(\d+)\s+unknowns\)"
                )
                m = re.search(percent_pattern, prop_line)
                if m:
                    current_record["percent"] = float(m.group(1))
                    current_record["align"] = int(m.group(2))
                    current_record["mismatch"] = int(m.group(3))
                    current_record["indel"] = int(m.group(4))
                    current_record["unknown"] = int(m.group(5))
                continue

    return result

def gmap_to_dataframe(data_dict):

    df = pd.DataFrame.from_dict(data_dict, orient='index')
    df = df.reset_index().rename(columns={'index': 'ref_name'})

    temp = df['ref_name'].str.split(r'\|', expand=True)

    df['gene_group'] = temp[0].astype(str) + "_group" + df['gene_number'].astype(str)
    df['ref_len'] = temp.iloc[:, -1].astype(int)

    df = df.drop('gene_number', axis=1)
    df = df.sort_values('roi_start')
    df = df[['ref_name', 'gene_group', 'ref_start', 'ref_end', 'roi_start', 'roi_end', 'strand', 'exons', 'percent', 'align', 'mismatch', 'indel','unknown','contig','ref_len']]
    
    return df

def main():
    if len(sys.argv) not in [2, 3]:
        print("Usage: python script.py input_file.txt [output_file.csv]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) == 3 else "output.csv"

    try:
        with open(input_file, 'r') as f:
            data = f.read()
    except Exception as e:
        print(f"Error reading file '{input_file}': {e}")
        sys.exit(1)

    parsed_data = parse_data(data)
    gmap_data = gmap_to_dataframe(parsed_data)
    
    gmap_data.to_csv(output_file, sep = '\t', index=False)
    print(f"Data written to {output_file}")

if __name__ == '__main__':
    main()
