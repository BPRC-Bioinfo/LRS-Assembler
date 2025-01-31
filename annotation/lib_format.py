import argparse
from Bio import SeqIO

# Argument parser setup
parser = argparse.ArgumentParser(description="Rename FASTA sequences by appending sequence length.")
parser.add_argument("input_fasta", help="Input FASTA file")
parser.add_argument("output_fasta", help="Output FASTA file")

# Parse the arguments
args = parser.parse_args()

# Open the input file and the output file for writing
with open(args.input_fasta, "r") as infile, open(args.output_fasta, "w") as outfile:
    # Iterate through each sequence in the input fasta file
    for record in SeqIO.parse(infile, "fasta"):
        # Get the sequence length
        seq_length = len(record.seq)
        # Create the new name by appending the length separated by a |
        new_id = f"{record.id}|{seq_length}"
        # Update the record ID
        record.id = new_id
        record.description = ""  # Optional: Clear the description to avoid duplication
        # Write the modified record to the output file
        SeqIO.write(record, outfile, "fasta")

print(f"Sequences have been renamed and saved to {args.output_fasta}")

