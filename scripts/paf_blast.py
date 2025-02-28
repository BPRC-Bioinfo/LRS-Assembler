# Blast genes from paf
# v0.0.1
# By Giang Le & Jaimy

import pandas as pd
import subprocess
import argparse
import shutil
import os
import re


def parse_line(line):
    try:
        ref_info, *rest = line.strip().split("\t")
        contig = rest[0]
        chr_start = int(rest[3])
        chr_end = int(rest[4])

        return ref_info, chr_start, chr_end, contig
    except ValueError:
        return None, None, None, None,

def gene_group(paf_file):
    groups = {}
    
    for line in paf_file:
        ref_info, chr_start, chr_end, contig = parse_line(line)

        if ref_info is not None:
            groups.setdefault(ref_info, []).append((chr_start, chr_end, contig))

    return groups

def blast_paf(groups, ref_file, fasta_file, threads, output_file):
    contigs_to_remove = set()
    for gene, values in groups.items():
        df = pd.DataFrame()
        for start, end, contig in values:
            contig_dir = os.path.join("map_temp", contig)
            contigs_to_remove.add(contig)
            os.makedirs(contig_dir, exist_ok=True)
            new_row = {'contig': contig, 'start': start, 'end': end, 'gene': gene}
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True) #efficiently append rows

        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)

        file_name = re.sub(r'[^a-zA-Z0-9]', '', gene.split('|')[0])
        bed_filename = os.path.join(contig_dir, f"{file_name}.bed")
        df.to_csv(bed_filename, sep='\t', header=False, index=False)

        chr_region_fasta = os.path.join(contig_dir, f"{file_name}_chr.fasta")
        try:        
            subprocess.run(["bedtools", "getfasta", "-fi", fasta_file, "-bed", bed_filename], stdout=open(chr_region_fasta, 'w'), check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running bedtools: {e}")
        except FileNotFoundError:
            print(f"Error: bedtools not found. Make sure it's installed and in your PATH.")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

        single_ref_fasta = os.path.join(contig_dir, f"{file_name}_ref.fasta")
        try:
            subprocess.run(["seqkit", "grep", "-p", gene, ref_file], stdout=open(single_ref_fasta, 'w'), check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running seqkit: {e}")
        except FileNotFoundError:
            print(f"Error: seqkit not found. Make sure it's installed and in your PATH.")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

        ref_log_file = os.path.join(contig_dir, f"{file_name}_db.log")
        try:    
            subprocess.run(["makeblastdb", "-in", single_ref_fasta, "-dbtype", "nucl", "-logfile", ref_log_file], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running makeblastdb: {e}")
        except FileNotFoundError:
            print(f"Error: makeblastdb not found. Make sure it's installed and in your PATH.")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

        try:
            with open(output_file, 'a') as outfile:
                subprocess.run(["blastn", "-query", chr_region_fasta, "-db", single_ref_fasta, "-num_threads", threads, "-word_size", "7", "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop qlen slen"], stdout=outfile, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running blastn: {e}")
        except FileNotFoundError:
            print(f"Error: blastn not found. Make sure it's installed and in your PATH.")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

    for contig in contigs_to_remove:  # Iterate through the *unique* contigs
        contig_dir = os.path.join("map_temp", contig)
        try:
            shutil.rmtree(contig_dir)
            print(f"Removed temporary directory: {contig_dir}")
        except FileNotFoundError:
            print(f"Temporary directory '{contig_dir}' not found (likely already removed).")
        except Exception as e:
            print(f"Error removing temporary directory {contig_dir}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract bed file from filtered and PAF files.",
                                     epilog="Example usage: python script_name.py -f filtered_file.tsv -p paf_file.paf -o output.bed")
    parser.add_argument("-f", "--fasta", type=str, help="Path to the FASTA input file.")
    parser.add_argument("-p", "--paf", type=str, help="Path to the PAF input file (TSV format).")
    parser.add_argument("-r", "--ref", type=str, help="Path to the reference input file (FASTA format).")
    parser.add_argument("-t", "--threads", default="4", type=str, help="Number of threads")
    parser.add_argument("-o", "--output", type=str, help="Path to the BED output file.")

    args = parser.parse_args()

    try:
        with open(args.paf, 'r') as f:
            paf_df = f.readlines()
    except FileNotFoundError:
        print("Error: One or all of the files specified were not found (even after initial existence check).")
        exit(1)
    except pd.errors.ParserError:
        print("Error: Could not parse one or both of the input files. Ensure they are valid TSV files.")
        exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading the input files: {e}")
        exit(1)

    groups = gene_group(paf_df)
    blast_paf(groups, args.ref, args.fasta, args.threads, args.output)



