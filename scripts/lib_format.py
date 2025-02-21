# Add sequence length to library
# v0.0.2
# By Giang Le

import argparse
from Bio import SeqIO

# Argument parser setup
parser = argparse.ArgumentParser(description="Rename FASTA sequences by appending sequence length.")
parser.add_argument("input_fasta", help="Input FASTA file")
parser.add_argument("output_fasta", help="Output FASTA file")

# Parse the arguments
args = parser.parse_args()

with open(args.input_fasta, "r") as infile, open(args.output_fasta, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        seq_length = len(record.seq)
        parts = record.id.split("|")

        if len(parts) > 1 and parts[-1].isdigit() and int(parts[-1]) == seq_length:
            SeqIO.write(record, outfile, "fasta")
            continue  # Skip to the next record

        new_id = f"{'|'.join(parts[:-1]) if len(parts) > 1 else record.id}|{seq_length}"
        record.id = new_id
        record.description = ""
        SeqIO.write(record, outfile, "fasta")

print(f"Library processed and saved to {args.output_fasta}")

