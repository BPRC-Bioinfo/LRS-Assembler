import argparse
import pandas as pd
import sys
import os
from os.path import commonprefix

def parse_args():
    parser = argparse.ArgumentParser(
        description='Load and combine cDNA, gDNA, and/or protein TSV files into one annotated TSV.'
    )
    parser.add_argument('--cdna', '-c', help='Path to the cDNA TSV file (or "None").')
    parser.add_argument('--gdna', '-g', help='Path to the gDNA TSV file (or "None").')
    parser.add_argument('--protein', '-p', help='Path to the protein TSV file (or "None").')
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Prefix for output file (default: common prefix of input file basenames or "combined")'
    )
    return parser.parse_args()

def load_tsv(path):
    try:
        return pd.read_csv(path, sep='\t', header=0, dtype=str)
    except Exception as e:
        sys.exit(f'Error: Could not read {path!r} as a TSV. Exception:\n  {e}')

def infer_output_prefix(paths):
    if not paths:
        return "combined"
    basenames = [os.path.splitext(os.path.basename(p))[0] for p in paths]
    cp = commonprefix(basenames).rstrip('_-')
    return cp or "combined"

def main():
    args = parse_args()

    # Gather and clean inputs
    provided = {
        'cdna': args.cdna,
        'gdna': args.gdna,
        'protein': args.protein
    }
    provided = {k: v for k, v in provided.items() if v and v != "None"}

    if not provided:
        sys.exit("Error: You must supply at least one of --cdna, --gdna, or --protein (and not 'None').")

    dfs = {}
    for kind, path in provided.items():
        if not os.path.isfile(path):
            sys.exit(f"Error: File not found: {path!r}")
        df = load_tsv(path)
        dfs[kind] = df

    # Combine into one DataFrame
    combined = pd.concat(dfs.values(), ignore_index=True)

    # Determine output path
    output_prefix = args.output if args.output and args.output != "None" else infer_output_prefix(provided.values())

    # Write output
    combined.to_csv(output_prefix, sep='\t', index=False)

    print(f"Combined {len(dfs)} files into: {output_prefix}")
    print(f"Rows written: {len(combined)}")

if __name__ == '__main__':
    main()
