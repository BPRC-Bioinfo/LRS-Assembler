# v0.0.7
# By Giang Le & Jaimy

import pandas as pd
import argparse
from dna_features_viewer import GraphicFeature, GraphicRecord
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import Range1d


def get_dna_viewer_plot(df_filtered: pd.DataFrame, output: str) -> dict:
    
    df_filtered = df_filtered.sort_values(1)
    if not df_filtered.empty:
        features = []
        for _, row in df_filtered.iterrows():
            color = "#ffcccc"  # Default color
            if "Flanking gene" in row[4]:
                color = "#ccffcc"  # Different color for flanking

            features.append(
                GraphicFeature(
                    start=int(row[1]),
                    end=int(row[2]),
                    strand=int(str(row[5]) + "1"),
                    label=f"{row[3]}",
                    color=color
                )
            )

        start_coord = df_filtered[1].min() - 100
        end_coord = df_filtered[2].max() - start_coord + 100
#        print (start_coord)
        record = GraphicRecord(
            sequence="ATCG",
            first_index=start_coord,
            sequence_length=end_coord,
            features=features
        )

        ax, _ = record.plot(figure_width=25)
       # ax.set_title(f"Annotation report for {sample_name}", loc='left', weight='bold')

        ax.figure.savefig(f'{output}.svg', bbox_inches='tight', dpi=600)
        ax.figure.savefig(f'{output}.pdf', bbox_inches='tight', dpi=600)
        
        zoom_start = int(df_filtered[1].min()) - 10000
        zoom_end = int(df_filtered[2].max()) + 10000

#        print (zoom_start, zoom_end)

        bokeh_plot = record.plot_with_bokeh(figure_width=20, figure_height=5)
        x_range = Range1d(zoom_start, zoom_end)

        bokeh_plot.x_range = x_range
        bokeh_plot.title.text_font_size = "20pt"
        bokeh_plot.title.align = "center"

        graph = file_html(bokeh_plot, CDN, "Region Viewer")

        return graph


def main():
    parser = argparse.ArgumentParser(description="Generate annotation flow plot from input file")
    parser.add_argument("-i", "--input", type=str, help="Path to the input Excel file")
    parser.add_argument("-o", "--output", type=str, help="Sample name for saving the figures")
    args = parser.parse_args()

    try:
        data = pd.read_csv(args.input, header = None, sep = "\t")
    except FileNotFoundError:
        print("Error: One or all of the files specified were not found (even after initial existence check).")
        exit(1)
    except pd.errors.ParserError:
        print("Error: Could not parse one or both of the input files. Ensure they are valid TSV files.")
        exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading the input files: {e}")
        exit(1)
    
    plot_html = get_dna_viewer_plot(data, args.output)



if __name__ == "__main__":
    main()



'''        
        features = [
            GraphicFeature(
                start=row["chr_start"],
                end=row["chr_end"],
                strand=int(str(row["strand"]) + "1"),
                label=f"{row['ref_name']}",
                color="#ffcccc"
            )
            for _, row in df_filtered.iterrows()
        ]
'''

