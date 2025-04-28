#!/usr/bin/env python3
import argparse
import pandas as pd

def process_flanks_intact(line, bed_file):
    fields = line.strip().split("\t")
    flank_genes, contig = fields[2:4]

    left_flank, right_flank = flank_genes.split("/")[0:2]  
    bed_file["flank"] = bed_file[3].str.split("_").str[-1]

    filtered = bed_file[(bed_file[0].astype(str) == contig) & (bed_file["flank"].isin([left_flank, right_flank]))]

    grouped_df = filtered.groupby("flank", as_index=False).agg({
    0: "first",
    1: "min",
    2: "max",
    3: "first",
    4: "first",
    5: "first"
    })

    return grouped_df

def process_flanks_fragmented(line, bed_file):
    fields = line.strip().split("\t")
    contig_list = fields[3].split('/')

    contigs = []
    for i in contig_list:
        contigs.append(i)

    return pd.DataFrame(contigs)



def main():
    parser = argparse.ArgumentParser(description="Process region intact bed coordinates.")
    parser.add_argument("-s", "--status", required=True, help="Input SAM file")
    parser.add_argument("-b", "--bed", required=True, help="Input BED file")
    parser.add_argument("-r", "--region", required=True, help="Region name")
    parser.add_argument("-o", "--output", help="Final output file")
    
    args = parser.parse_args()

    bed_df = pd.read_csv(args.bed, sep="\t", header=None)

    final_grouped_df = []
    with open(args.status) as fstatus:
        for line in fstatus:
            if "closed" in line and "flanks" in line:
                grouped_df = process_flanks_intact(line, bed_df)
            elif "fragmented" in line:
                grouped_df = process_flanks_fragmented(line, bed_df)
            
            final_grouped_df.append(grouped_df)

    if final_grouped_df:
        final_flanks = pd.concat(final_grouped_df, ignore_index=True)
        if final_flanks.shape[1] == 7:
            final_bed = final_flanks.groupby(0, as_index=False).agg({
            0: "first",
            1: "min",
            2: "max",
            3: "first",
            4: "first",
            5: "first"
            })

            final_bed[3] = args.region
            # Condition for length
            # Add 50 bp in front and end
    #        final_bed[1] = final_bed[1] - 50
    #        final_bed[2] = final_bed[2] + 50
            final_bed.to_csv(f"{args.output}", sep="\t", index=False, header=False)
        else:
            final_flanks.to_csv(f"{args.output}", sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()
