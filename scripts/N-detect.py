import sys

def detect_N_locations(fasta_file):
    with open(fasta_file, 'r') as f:
        contig_name = ""
        contig_sequence = ""
        for line in f:
            if line.startswith(">"):
                # If we have previously read a contig, process it
                if contig_name and contig_sequence:
                    detect_N(contig_name, contig_sequence)
                    contig_name = ""
                    contig_sequence = ""
                contig_name = line.strip()[1:]
            else:
                contig_sequence += line.strip()
        # Process the last contig in the file
        if contig_name and contig_sequence:
            detect_N(contig_name, contig_sequence)

def detect_N(contig_name, sequence):
    n_positions = []
    start = None
    end = None
    length = len(sequence)
    for i, base in enumerate(sequence):
        if base == 'N':
            if start is None:
                start = i
            end = i
        elif start is not None:
            n_positions.append((start, end, end - start + 1))
            start = None
            end = None
    if n_positions:
        for pos in n_positions:
            print(f"{contig_name} Start: {pos[0]}, End: {pos[1]}, Length: {pos[2]}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_fasta_file>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    detect_N_locations(fasta_file)

