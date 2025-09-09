#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO

def load_map(map_file):
    df = pd.read_csv(map_file, sep=None, engine='python')
    if 'accession' not in df.columns or 'chromosome' not in df.columns:
        raise ValueError("Mapping file must contain 'accession' and 'chromosome' columns")
    return df.set_index('accession')['chromosome'].to_dict()

def rename_fasta(input_fasta, output_fasta, acc2chr, pattern):
    total = renamed = 0

    with open(output_fasta, 'w') as out_handle:
        for rec in SeqIO.parse(input_fasta, 'fasta'):
            total += 1
            orig_id = rec.id
            acc = orig_id.split()[0]

            if acc in acc2chr:
                chr_name = acc2chr[acc]
                # New header is: <pattern>_chr<chromosome>
                rec.id = f"{pattern}_chr{chr_name}"
                rec.description = ""  # clear description for cleanliness
                renamed += 1
            else:
                # Keep original header exactly as-is
                rec.id = orig_id
                # rec.description left untouched if more details are desired

            SeqIO.write(rec, out_handle, 'fasta')

    print(f"# Done. {renamed}/{total} sequences renamed; all {total} sequences retained.")

def main():
    parser = argparse.ArgumentParser(description="Rename FASTA headers using accession mapping.")
    parser.add_argument("--ac_file", required=True, help="CSV/TSV with 'accession' and 'chromosome' columns")
    parser.add_argument("--input_fasta", required=True, help="Original FASTA input")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA")
    parser.add_argument(
        "--header_pattern", required=True,
        help="New header pattern, e.g. '{sample}' or '{sample}_hap1' "
    )

    args = parser.parse_args()
    acc2chr = load_map(args.ac_file)
    rename_fasta(args.input_fasta, args.output_fasta, acc2chr, args.header_pattern)

if __name__ == "__main__":
    main()
