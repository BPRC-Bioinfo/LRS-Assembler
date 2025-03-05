# Add sequence length to library
# v0.0.3
# By Giang Le

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import sys


parser = argparse.ArgumentParser(description="Rename FASTA sequences by appending sequence length.")
parser.add_argument("input_fasta", help="Input FASTA file")

args = parser.parse_args()

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

'''
def remove_dashes(sequence):
  seq_obj = Seq(sequence)
  cleaned_seq = str(seq_obj).replace("-", "")
  return cleaned_seq

def process_fasta(input_fasta):

    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            parts = record.id.split("|")

            # Remove dashes from the sequence
            cleaned_sequence = remove_dashes(str(record.seq))
            record.seq = Seq(cleaned_sequence) #update the seq object
            seq_length = len(record.seq) #update the length after dash removal

            if len(parts) > 1 and parts[-1].isdigit() and int(parts[-1]) == seq_length:
                SeqIO.write(record, sys.stdout, "fasta")
                continue  # Skip to the next record

            new_id = f"{'|'.join(parts[:-1]) if len(parts) > 1 else record.id}|{seq_length}"
            record.id = new_id
            record.description = ""

            SeqIO.write(record, sys.stdout, "fasta")
'''
