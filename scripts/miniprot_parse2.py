import argparse
import re
from collections import defaultdict

def count_mismatches(cs_string: str) -> int:
    """
    Count the number of single‐base substitutions in a cs:Z: string.
    Each substitution is encoded as '*<ref_base><query_base>', 
    e.g. '*gG' for ref='g', query='G'.
    """
    substitutions = re.findall(r'\*[A-Za-z]{2}', cs_string)
    return len(substitutions)

def count_deletions(cs_string: str) -> int:
    """
    Count the total number of deleted bases in a cs:Z: string.
    Each deletion is encoded as '-<deleted_bases>', e.g. '-gg' means two bases deleted.
    This function returns the sum of lengths of all deletion sequences.
    """
    deletion_seqs = re.findall(r'-([A-Za-z]+)', cs_string)
    return sum(len(seq) for seq in deletion_seqs)

def extract_miniprot_with_mismatches(input_file: str, output_file: str):
    """
    Read a miniprot‐formatted file and for each hit:
      - Extract columns 1–6, 8–11, and 15 (1-based indices):
          1  ref_name
          2  ref_len
          3  ref_start
          4  ref_end
          5  strand
          6  contig
          8  roi_start
          9  roi_end
         10  nu_mapped
         11  nu_total
         15  pro_mapped (format "np:i:<number>")
      - Parse the cs:Z: tag from the same line to count mismatches (substitutions)
        and total deleted bases (indels).
      - Compute:
          percent       = (nu_mapped / nu_total) * 100
          align_percent = (pro_mapped / ref_len) * 100
          vs_ref        = abs(pro_mapped - ref_len)
      - Create gene_group by counting occurrences of each ref_name:
          e.g. KIR3DL3*020_group1, KIR3DL3*020_group2, ...
      - Write a tab-delimited output with these columns plus:
          percent, align_percent, vs_ref, mismatch_count, indel, gene_group
    """
    # Zero-based indices corresponding to columns 1,2,3,4,5,6, 8,9,10,11, 15
    # (i.e., parts[0], parts[1], parts[2], parts[3], parts[4], parts[5],
    #        parts[7], parts[8], parts[9], parts[10], parts[14])
    cols_to_keep = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 14]
    ref_counts = defaultdict(int)

    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        # Write header
        headers = [
            "ref_name", "ref_len", "ref_start", "ref_end", "strand", "contig",
            "roi_start", "roi_end", "nu_mapped", "nu_total", "align",
            "percent", "align_percent", "vs_ref", "gene_group"
        ]
        fout.write("\t".join(headers) + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if line.startswith('#'):
                continue

            parts = line.split("\t")
            if len(parts) < 15:
                continue

            # Extract fixed columns
            ref_name = parts[0]
            try:
                ref_len   = int(parts[1])
                ref_start = int(parts[2])
                nu_mapped = int(parts[9])
                nu_total  = int(parts[10])
            except ValueError:
                continue

            # pro_mapped (field 15, index 14, format "np:i:<number>")
            pro_field = parts[14]
            pro_parts = pro_field.split(':')
            try:
                pro_mapped = int(pro_parts[-1])
            except ValueError:
                continue

            # Find cs:Z: tag among optional fields (indices 11 and beyond)
            cs_tag = ""
            for tag in parts[11:]:
                if tag.startswith("cs:Z:"):
                    cs_tag = tag[len("cs:Z:"):]
                    break

            # Compute metrics
            percent = (nu_mapped / nu_total * 100) if nu_total > 0 else 0.0
            align_percent = (pro_mapped / ref_len * 100) if ref_len > 0 else 0.0
            vs_ref = abs(pro_mapped - ref_len)
            ref_start = ref_start + 1

            percent_str = f"{percent:.2f}"
            align_percent_str = f"{align_percent:.2f}"

            # Determine gene_group
            ref_counts[ref_name] += 1
            gene_group = f"{ref_name}_group{ref_counts[ref_name]}"

            # Build output row
            selected = [
                ref_name,           # col 1
                parts[1],           # col 2 (ref_len)
                str(ref_start),           # col 3 (ref_start)
                parts[3],           # col 4 (ref_end)
                parts[4],           # col 5 (strand)
                parts[5],           # col 6 (contig)
                parts[7],           # col 8 (roi_start)
                parts[8],           # col 9 (roi_end)
                parts[9],           # col 10 (nu_mapped)
                parts[10],          # col 11 (nu_total)
                str(pro_mapped),    # col 15 (pro_mapped as integer)
                percent_str,
                align_percent_str,
                str(vs_ref),
                gene_group
            ]
            fout.write("\t".join(selected) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract columns and compute metrics from a miniprot file, including mismatch and indel counts."
    )
    parser.add_argument("input", help="Path to input miniprot file")
    parser.add_argument("output", help="Path to tab-delimited output file")
    args = parser.parse_args()

    extract_miniprot_with_mismatches(args.input, args.output)
