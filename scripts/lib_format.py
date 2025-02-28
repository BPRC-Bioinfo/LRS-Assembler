# Add sequence length to library
# v0.0.3
# By Giang Le

import argparse
from Bio import SeqIO
import sys

parser = argparse.ArgumentParser(description="Rename FASTA sequences by appending sequence length.")
parser.add_argument("input_fasta", help="Input FASTA file")

args = parser.parse_args()

with open(args.input_fasta, "r") as infile:
    for record in SeqIO.parse(infile, "fasta"):
        seq_length = len(record.seq)
        parts = record.id.split("|")

        if len(parts) > 1 and parts[-1].isdigit() and int(parts[-1]) == seq_length:

            SeqIO.write(record, sys.stdout, "fasta")
            continue  # Skip to the next record

        new_id = f"{'|'.join(parts[:-1]) if len(parts) > 1 else record.id}|{seq_length}"
        record.id = new_id
        record.description = ""

        SeqIO.write(record, sys.stdout, "fasta")
