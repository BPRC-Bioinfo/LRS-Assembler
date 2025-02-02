# ver 0.0.9
# By: Giang Le & Jaimy
# Process paf with blastn information 

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

    # Process the number of subgroups
    result_data = []
    dict_data = {}

    # Now process each group
    for ref_name, group_lines in groups.items():
        # Processing multiple hits
#        print (group_lines)
#        print ("number of small list")
#        print (len(group_lines))
        current_ref_start = None
        subgroup = 0
        dict_data[ref_name]= {}
        for i, (ref_start_coord, ref_end_coord, chr_start_coord, chr_end_coord, length, original_line) in enumerate(group_lines):
            if current_ref_start is None:  # First line of the group
                current_ref_start = ref_start_coord
                subgroup += 1
#                print ("First current line and subgroup1")
#                print (subgroup)
#                print (original_line)
                dict_data[ref_name][subgroup]= [original_line]

            # Processing group with multiple hits
            elif ref_start_coord == current_ref_start or ref_start_coord < current_ref_start:
                subgroup += 1
#                print ("Start of new subgroup")
#                print (subgroup)
#                print (original_line)
                dict_data[ref_name][subgroup]= [original_line]
            else:
#                print ("Add lines to current group")
#                print (subgroup)
#                print (original_line)
#                if subgroup not in dict_data[ref_name]:
#                    print ("When does this trigger?")
#                    dict_data[ref_name][subgroup] = []  # Initialize list if subgroup doesn't exist
                dict_data[ref_name][subgroup].append(original_line)

            # Update current_ref_start for the next iteration
            current_ref_start = ref_start_coord
    
#    print(dict_data)

    result_data = []

    # Loop through the outer dictionary (samples)
    for ref_name, subgroups in dict_data.items():
        # Loop through the nested dictionary (groups)
        for subgroup_name, subgroup_info in subgroups.items():
#            print ([item.strip("\n").split(" ") for item in subgroup_info])
            ref_min_start = min(item.strip("\n").split(" ")[1] for item in subgroup_info)
            ref_max_end = max(item.strip("\n").split(" ")[2] for item in subgroup_info)
            chr_min_start = min(int(item.strip("\n").split(" ")[3]) for item in subgroup_info)
            chr_max_end = max(int(item.strip("\n").split(" ")[4]) for item in subgroup_info)
            chr_identity = mean(float(item.strip("\n").split(" ")[11]) 
                    if item.strip("\n").split(" ")[11].replace('.', '', 1).isdigit() 
                    else 0
                    for item in subgroup_info)
            chr_coverage = sum(int(item.strip("\n").split(" ")[12]) 
                    if item.strip("\n").split(" ")[12].replace('.', '', 1).isdigit() 
                    else 0
                    for item in subgroup_info)

            chr_coverage_percent = chr_coverage / int(ref_name.split("|")[1]) * 100
            missing_fragments = sum(1 if not item.strip("\n").split(" ")[11].strip() else 0 for item in subgroup_info)
            chr_insertion = sum(int(item.strip("\n").split(" ")[8]) for item in subgroup_info)
            chr_deletion = sum(int(item.strip("\n").split(" ")[7]) for item in subgroup_info)
            chr_mismatch = sum(int(item.strip("\n").split(" ")[13]) 
                   if item.strip("\n").split(" ")[13].isdigit() 
                   else 0
                   for item in subgroup_info)
            chr_gap = sum(int(item.strip("\n").split(" ")[14]) 
                   if item.strip("\n").split(" ")[14].isdigit() 
                   else 0
                   for item in subgroup_info)
            strand = [item.strip("\n").split(" ")[9] for item in subgroup_info]
            if all(s == strand[0] for s in strand):
                common_strand = strand[0]  # If they are all the same, pick the first one
            else:
                common_strand = None

            result_data.append({
                'Group_name':ref_name,
                'ref_start': ref_min_start,
                'ref_end': ref_max_end,
                'chr_start': chr_min_start,
                'chr_end': chr_max_end,
                'identity_percent': f"{chr_identity:.2f}",
                'aligned_percent': f"{chr_coverage_percent:.2f}",
                'miss_exons': missing_fragments,
                'minimap_ins': chr_insertion,
                'minimap_del': chr_deletion,
                'blast_mismatch': chr_mismatch,
                'blast_gap': chr_gap,
                'strand': common_strand
            })

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