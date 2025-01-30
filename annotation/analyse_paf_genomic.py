# ver 0.0.8
# By: Giang Le & Jaimy

import argparse
import pandas as pd
from statistics import mean


def parse_line(line):
    parts = line.strip().split(" ")
    ref_name = parts[0]
    ref_start_coord = int(parts[1])
    ref_end_coord = int(parts[2])
    chr_start_coord = int(parts[3])
    chr_end_coord = int(parts[4])
    details = parts[0].split("|")
    length = int(details[-1])

    return ref_name, ref_start_coord, ref_end_coord, chr_start_coord, chr_end_coord, length, line


def group_by_ref_name(file_path):
    groups = {}

    # Open the file and read lines
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Iterate over each line and group by the first column
    for line in lines:
        ref_name, ref_start_coord, ref_end_coord, chr_start_coord, chr_end_coord, length, original_line = parse_line(line)
        
        # Create a dictionary of groups by the first column
        if ref_name not in groups:
            groups[ref_name] = []
        groups[ref_name].append((ref_start_coord, ref_end_coord, chr_start_coord, chr_end_coord, length, original_line))
#    print (groups)


    # Initialize the list for the dataframe data
    result_data = []

    # Now process each group
    for ref_name, group_lines in groups.items():
        # Initialize subgroups
        subgroups = []
        current_subgroup = []

        # Iterate over the lines and create subgroups while preserving the original order
        for i, (ref_start_coord, ref_end_coord, chr_start_coord, chr_end_coord, length, original_line) in enumerate(group_lines):
            # If the next ref_start_coord is smaller, start a new subgroup
            if i > 0 and ref_start_coord < group_lines[i - 1][0]:
                subgroups.append(current_subgroup)
                current_subgroup = []
            current_subgroup.append((ref_start_coord, ref_end_coord, chr_start_coord, chr_end_coord, length, original_line))

        # Append the last subgroup
        if current_subgroup:
            subgroups.append(current_subgroup)

#        print (ref_name)
#        print (subgroups)

        # Process subgroups and calculate aligned length
        for idx, subgroup in enumerate(subgroups):

#            for sublist in subgroup:
#                entry_count = len(sublist)
#                print (entry_count) 

            # Calculate alignment length if BLAST results present
            aligned_length = sum(int(parts[6]) if len((parts := last.strip("\n").split(" "))) > 10 else 0 for _, _, _, _, _, last in subgroup) / length *100
#            print (aligned_length)

            # Skip processing the whole subgroup if aligned_length exceeds 100
            if aligned_length > 100:
                # Loop through each line in the subgroup
                for line_idx, (line_ref_start_coord, line_ref_end_coord, line_chr_start_coord, line_chr_end_coord, _, line_last) in enumerate(subgroup):
                    line_parts = line_last.strip("\n").split(" ")

                    # Calculate individual alignment length for each line
                    line_aligned_length = int(line_parts[6]) / length * 100 if len(line_parts) > 10 else 0

                    # Only add lines with aligned_percent <= 100
                    if line_aligned_length <= 100:
                        # Calculate identity, mismatches, gaps, etc.
                        line_identity = float(line_parts[10]) if len(line_parts) > 10 else 0
                        line_mismatch = int(line_parts[12]) if len(line_parts) > 10 else int(line_parts[2]) - int(line_parts[1])
                        line_gap = int(line_parts[13]) if len(line_parts) > 10 else int(line_parts[2]) - int(line_parts[1])
                        line_insertion = int(line_parts[8]) if len(line_parts) > 10 else 0
                        line_deletion = int(line_parts[7]) if len(line_parts) > 10 else 0
                        line_missing_fragment = 1 if len(line_parts) <= 10 else 0

                        # Create a unique name for each line and append it to the result
                        line_ref_name = ref_name.split("|")[0]
                        result_data.append({
                            'Group_name': f"{line_ref_name}_like{line_idx+1}|{length}",
                            'Subgroup': line_idx + 1,
                            'ref_start': line_ref_start_coord,
                            'ref_end': line_ref_end_coord,
                            'chr_start': line_chr_start_coord,
                            'chr_end': line_chr_end_coord,
                            'identity_percent': f"{line_identity:.2f}",
                            'aligned_percent': f"{line_aligned_length:.2f}",
                            'missing_exons': line_missing_fragment,
                            'mismatch_sum': line_mismatch,
                            'gap_sum': line_gap,
                            'insertions': line_insertion,
                            'deletion': line_deletion,
                        })

            else:  # Process whole subgroup if aligned_length <= 100
                min_ref_start_coord = min(ref_start_coord for ref_start_coord, _, _, _, _, _ in subgroup)
                max_ref_end_coord = max(ref_end_coord for _, ref_end_coord, _, _, _, _ in subgroup)
                min_chr_start_coord = min(chr_start_coord for _, _, chr_start_coord, _, _, _ in subgroup)
                max_chr_end_coord = max(chr_end_coord for _, _, _, chr_end_coord, _, _ in subgroup)

                identity_avg = mean(float(parts[10]) if len((parts := last.strip("\n").split(" "))) > 10 else 0 for _, _, _, _, _, last in subgroup)

                missing_fragments = sum(1 for _, _, _, _, _, last in subgroup if len((parts := last.strip("\n").split(" "))) <= 10)
                mismatch_sum = sum(int(parts[12]) if len((parts := last.strip("\n").split(" "))) > 10 else int(parts[2]) - int(parts[1]) for _, _, _, _, _, last in subgroup)
                gap_sum = sum(int(parts[13]) if len((parts := last.strip("\n").split(" "))) > 10 else int(parts[2]) - int(parts[1]) for _, _, _, _, _, last in subgroup)
                insertion_sum = sum(int(last.strip("\n").split(" ")[8]) for _, _, _, _, _, last in subgroup)
                deletion_sum = sum(int(last.strip("\n").split(" ")[7]) for _, _, _, _, _, last in subgroup)

                result_data.append({
                    'Group_name': ref_name,
                    'Subgroup': idx + 1,
                    'ref_start': min_ref_start_coord,
                    'ref_end': max_ref_end_coord,
                    'chr_start': min_chr_start_coord,
                    'chr_end': max_chr_end_coord,
                    'identity_percent': f"{identity_avg:.2f}",
                    'aligned_percent': f"{aligned_length:.2f}",
                    'missing_exons': missing_fragments,
                    'mismatch_sum': mismatch_sum,
                    'gap_sum': gap_sum,
                    'insertions': insertion_sum,
                    'deletion': deletion_sum,
                })


    # Create the dataframe
    df = pd.DataFrame(result_data)
    df = df.sort_values('chr_start')

    return df

if __name__ == "__main__":
    # Set up argparse to handle command-line arguments
    parser = argparse.ArgumentParser(description="Group data by the first column and split into sub-groups based on fragment length and coordinates.")
    parser.add_argument("input_file", help="Input file containing the data to process")
    parser.add_argument("output_file", help="Output file name (without extension)")

    # Parse arguments
    args = parser.parse_args()

    # Call the group_by_ref_name function with the input fileqqq
    df = group_by_ref_name(args.input_file)

    # Save the dataframe to CSV and Excel files with the specified output name
    df.to_csv(f"{args.output_file}.csv", sep='\t', index=False) 
    df.to_excel(f"{args.output_file}.xlsx", index=False, sheet_name="Full_info")

    # Display the filtered dataframe
    print(df.to_string())