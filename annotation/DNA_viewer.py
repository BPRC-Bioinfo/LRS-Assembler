import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
import argparse


def plot(df, sample_name) -> None:
    features = [
        GraphicFeature(
            start=row["chr_start"],
            end=row["chr_end"],
            label=f"{row['Simple name']}"
        )
        for _, row in df.iterrows()
    ]

    start_coord = df["chr_start"].min() - 10000
    end_coord = df["chr_end"].max() - start_coord + 10000

    record = GraphicRecord(first_index=start_coord, sequence_length=end_coord, features=features)

    ax, _ = record.plot(figure_width=25)
    ax.set_title(f"Annotation report for {sample_name}", loc='left', weight='bold')

    # Save figures with the sample name included in the filenames
    ax.figure.savefig(f'figure/annotation_flow_{sample_name}.png', bbox_inches='tight', dpi=600)
    ax.figure.savefig(f'figure/annotation_flow_{sample_name}.svg', bbox_inches='tight', dpi=600)
    ax.figure.savefig(f'figure/annotation_flow_{sample_name}.pdf', bbox_inches='tight', dpi=600)


def main():
    parser = argparse.ArgumentParser(description="Generate annotation flow plot from input file")
    parser.add_argument("input_file", type=str, help="Path to the input Excel file")
    parser.add_argument("sample_name", type=str, help="Sample name for saving the figures")
    args = parser.parse_args()

    # Read the Excel file using the provided input path
    df = pd.read_excel(args.input_file)

    # Process the dataframe as before
    df["Simple name"] = df["Group_name"].str.split("|").str[0]

    # Call the plot function with the sample name
    plot(df, args.sample_name)


if __name__ == "__main__":
    main()

