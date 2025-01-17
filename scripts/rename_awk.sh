#!/bin/bash

# Ensure that the script receives two arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <keep_list_file> <sample_name> <input_fasta_file> <output_fasta_file>"
    exit 1
fi

# Input files
KEEP_LIST_FILE="$1"
SAMPLE_NAME="$2"
INPUT_FASTA="$3"
OUTPUT_FASTA="$4"

# Read the keep list file and generate the AWK lookup pattern
KEEP_PATTERN=$(awk '{printf "|%s", $0}' "$KEEP_LIST_FILE")
KEEP_PATTERN=${KEEP_PATTERN:1} # Remove leading '|'

# Counter for renaming other contigs
counter=1

# AWK command to rename contigs
awk -v sample="$SAMPLE_NAME" -v keep_pattern="$KEEP_PATTERN" -v counter="$counter" '
  BEGIN { 
    # Convert the keep pattern into an array for easy lookup
    n = split(keep_pattern, keep_names, "|");
    for (i = 1; i <= n; i++) {
      keep[keep_names[i]] = 1;
    }
  }
  /^>/ {
    # Check if the header (contig name) is in the keep list
    header_name = substr($0, 2);
    if (header_name in keep) {
      # Print the original header
      print $0;
    } else {
      # Print the new header
      print ">"sample"_contig_" counter;
      counter++;
    }
  }
  !/^>/ {
    # Print the sequence lines as they are
    print $0;
  }
' "$INPUT_FASTA" > "$OUTPUT_FASTA"

echo "Output written to $OUTPUT_FASTA"

