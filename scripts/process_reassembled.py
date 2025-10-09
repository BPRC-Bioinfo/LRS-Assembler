#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess

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
    parser.add_argument("-b", "--bed", required=True, help="Path to the BED file for the flanking genes")
    parser.add_argument("-f", "--fasta", required=True, help="Path to reference genome FASTA file")
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

    gene_coords = {}  # (gene, contig) -> (start, end)

    with open(args.bed, "r") as floc:
        for line in floc:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            contig = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene_tag = parts[3]
            gene_name = gene_tag.split("_")[-1]

            if left_flank and gene_name == left_flank:
                left_contigs.append(contig)
                gene_coords[(left_flank, contig)] = (start, end)
            elif right_flank and gene_name == right_flank:
                right_contigs.append(contig)
                gene_coords[(right_flank, contig)] = (start, end)

    left_contigs = sorted(set(left_contigs))
    right_contigs = sorted(set(right_contigs))

    output_lines = []
    left_set = set(left_contigs)
    right_set = set(right_contigs)

    closed_regions = {}         # contig -> (start-1, end)
    fragmented_contigs = set()  # unique contigs from fragmented pairs (no *)

    if left_flank and right_flank:
        common = left_set & right_set
        if common:
            for contig in sorted(common):
                output_lines.append(f"flanks\tclosed\t{left_flank}/{right_flank}\t{contig}")
                left_pos = gene_coords.get((left_flank, contig))
                right_pos = gene_coords.get((right_flank, contig))
                if left_pos and right_pos:
                    region_start = min(left_pos[0], right_pos[0]) - 1
                    region_end = max(left_pos[1], right_pos[1])
                    closed_regions[contig] = (region_start, region_end)

        # Fragmented
        fragmented_left = left_set - common
        fragmented_right = right_set - common
        if fragmented_left or fragmented_right:
            left_str = " ".join(sorted(fragmented_left)) if fragmented_left else "*"
            right_str = " ".join(sorted(fragmented_right)) if fragmented_right else "*"
            output_lines.append(f"flanks\tfragmented\t{left_flank}/{right_flank}\t{left_str}/{right_str}")

            if fragmented_left and fragmented_right:
                for l in sorted(fragmented_left):
                    for r in sorted(fragmented_right):
                        fragmented_contigs.update([l, r])
            elif fragmented_left:
                fragmented_contigs.update(fragmented_left)
            elif fragmented_right:
                fragmented_contigs.update(fragmented_right)

        if not left_set and not right_set:
            output_lines.append("flanks\tmissing\t*\t*")

    elif left_flank:
        if left_set:
            for contig in sorted(left_set):
                output_lines.append(f"telomere\tclosed\t{left_flank}/\t{contig}")
        else:
            output_lines.append(f"telomere\tmissing\t{left_flank}/\t*")

    elif right_flank:
        if right_set:
            for contig in sorted(right_set):
                output_lines.append(f"telomere\tclosed\t/{right_flank}\t{contig}")
        else:
            output_lines.append(f"telomere\tmissing\t/{right_flank}\t*")

    else:
        output_lines.append("flanks\tmissing\t*\t*")

    # Write main status file
    unique_lines = remove_duplicates(output_lines)

    roi_reassembling_status = args.out + "_missing_status.txt"
    with open(roi_reassembling_status, "w") as fstatus:
        for line in unique_lines:
            fstatus.write(line + "\n")

    # ---- Output closed regions as BED ----
    bed_file = args.out + ".closed.bed"
    with open(bed_file, "w") as b:
        for contig, (start, end) in closed_regions.items():
            b.write(f"{contig}\t{start}\t{end}\n")

    # ---- Output fragmented contigs for seqkit ----
    frag_file = args.out + ".fragmented.txt"
    with open(frag_file, "w") as f:
        for contig in sorted(fragmented_contigs):
            f.write(contig + "\n")

    # --- Step 1: Extract closed regions using bedtools ---
    combined_fasta = args.out + ".fa"
    bedtools_cmd = ["bedtools", "getfasta", "-fi", args.fasta, "-bed", bed_file, "-fo", combined_fasta]
    print("\nRunning:", " ".join(bedtools_cmd))
    subprocess.run(bedtools_cmd, check=True)

    # --- Step 2: Extract fragmented contigs using seqkit, append to combined FASTA ---
    seqkit_cmd = ["seqkit", "grep", "-f", frag_file, args.fasta]
    print("Running (and appending to FASTA):", " ".join(seqkit_cmd))

    with open(combined_fasta, "a") as fout:
        proc = subprocess.run(seqkit_cmd, check=True, stdout=fout)

    print("\nCombined closed + fragmented FASTA written to:", combined_fasta)

    try:
        os.remove(bed_file)
        os.remove(frag_file)
        print("Cleaned up temporary files:", bed_file, frag_file)
    except Exception as e:
        print("Warning: Could not remove temp files:", str(e))

if __name__ == "__main__":
    main()
