#!/usr/bin/env python3
import sys
import argparse

def remove_duplicates(lines):
    seen = set()
    unique = []
    for line in lines:
        if line not in seen:
            seen.add(line)
            unique.append(line)
    return unique

def main():
    parser = argparse.ArgumentParser(description="Check flanking genes")
    parser.add_argument("-i", "--info", required=True, help="Path to the info file")
    parser.add_argument("-l", "--loc", required=True, help="Path to the flanking genes loc file")
    parser.add_argument("-o", "--out", required=True, help="Path to final status file")

    args = parser.parse_args()

    left_flank = ""
    right_flank = ""
    with open(args.info, "r") as finfo:
        for line in finfo:
            line = line.strip()
            if line.startswith("left:"):
                left_flank = line.split(":", 1)[1].strip()
            elif line.startswith("right:"):
                right_flank = line.split(":", 1)[1].strip()

    left_contigs = []
    right_contigs = []

    with open(args.loc, "r") as floc:
        for line in floc:
            line = line.strip()
            if not line:
                continue
            # Expected tab-separated fields; we use at least the first and third.
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            region_info = parts[0]
            contig = parts[2]
            # Extract gene_name: last element after splitting on underscore.
            gene_name = region_info.split("_")[-1]
            if left_flank and gene_name == left_flank:
                left_contigs.append(contig)
            elif right_flank and gene_name == right_flank:
                right_contigs.append(contig)

    left_contigs = sorted(set(left_contigs))
    right_contigs = sorted(set(right_contigs))
    
    output_lines = []

    left_set = set(left_contigs)
    right_set = set(right_contigs)

    if left_flank and right_flank:
        common = left_set & right_set
        if common:
            for contig in sorted(common):
                output_lines.append(f"flanks\tclosed\t{left_flank}/{right_flank}\t{contig}")
        elif not left_set and not right_set:
            output_lines.append("flanks\tmissing\t*\t*")
        else:
            left_str = " ".join(sorted(left_set)) if left_set else "*"
            right_str = " ".join(sorted(right_set)) if right_set else "*"
            output_lines.append(f"flanks\tfragmented\t{left_flank}/{right_flank}\t{left_str}/{right_str}")
    elif left_flank:
        for contig in sorted(left_set):
            output_lines.append(f"telomere\tclosed\t{left_flank}/\t{contig}")
    elif right_flank:
        for contig in sorted(right_set):
            output_lines.append(f"telomere\tclosed\t/{right_flank}\t{contig}")

    # Remove duplicates (if any) and write final status file
    unique_lines = remove_duplicates(output_lines)

    with open(args.out, "w") as fstatus:
        for line in unique_lines:
            fstatus.write(line + "\n")

if __name__ == "__main__":
    main()
