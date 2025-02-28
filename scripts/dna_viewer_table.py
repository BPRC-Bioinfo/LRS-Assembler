import pandas as pd
import argparse
from dna_features_viewer import GraphicFeature, GraphicRecord
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import Range1d

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)



def get_dna_viewer_plot(df_filtered: pd.DataFrame, output: str) -> dict:
    if not df_filtered.empty:
        features = []
        for _, row in df_filtered.iterrows():
            color = "#ffcccc"  # Default color
            if "Flanking gene" in row["blast_percent"]:
                color = "#ccffcc"  # Different color for flanking

            features.append(
                GraphicFeature(
                    start=row["chr_start"],
                    end=row["chr_end"],
                    strand=int(str(row["strand"]) + "1"),
                    label=f"{row['ref_name']}",
                    color=color
                )
            )

        start_coord = df_filtered["chr_start"].min() - 100
        end_coord = df_filtered["chr_end"].max() - start_coord + 100

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
        
        zoom_start = int(df_filtered['chr_start'].min()) - 10000
        zoom_end = int(df_filtered['chr_end'].max()) + 10000

        print (zoom_start, zoom_end)

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
        data = pd.read_csv(args.input, sep = "\t")
    except FileNotFoundError:
        print("Error: One or all of the files specified were not found (even after initial existence check).")
        exit(1)
    except pd.errors.ParserError:
        print("Error: Could not parse one or both of the input files. Ensure they are valid TSV files.")
        exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading the input files: {e}")
        exit(1)
    # Read the Excel file using the provided input path
    
    plot_html = get_dna_viewer_plot(data, args.output)

    table_html = f"""
    <div class="table-responsive" style="height: 750px; overflow-y: scroll;">
        <table class="table table-hover mt-3" id="dataTable">
            <thead class="table-light">
                <tr>
                    {''.join(f'<th>{col}</th>' for col in data.columns if col in ['ref_name', 'chr_start', 'chr_end', 'blast_percent', 'blast_align', 'blast_mismatch', 'blast_gap', 'strand', 'ref_len'])}
                </tr>
            </thead>
            <tbody>
                {''.join('<tr>' + ''.join(f'<td>{row[col]}</td>' for col in data.columns if col in ['ref_name', 'chr_start', 'chr_end', 'blast_percent', 'blast_align', 'blast_mismatch', 'blast_gap', 'strand', 'ref_len'])
                         + '</tr>' for _, row in data.iterrows())}
            </tbody>
        </table>
    </div>
    """


    html_template = f"""
    <!DOCTYPE html>
    <html lang="nl">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Report</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    </head>
     <script>
            function filterTable() {{
                var input, filter, table, tr, td, i, j, txtValue, count;
                input = document.getElementById("searchInput");
                filter = input.value.toUpperCase();
                table = document.getElementById("dataTable");
                tr = table.getElementsByTagName("tr");
                count = 0;

                for (i = 1; i < tr.length; i++) {{
                    tr[i].style.display = "none";
                    td = tr[i].getElementsByTagName("td");
                    for (j = 0; j < td.length; j++) {{
                        if (td[j]) {{
                            txtValue = td[j].textContent || td[j].innerText;
                            if (txtValue.toUpperCase().indexOf(filter) > -1) {{
                                tr[i].style.display = "";
                                count++;
                                break;
                            }}
                        }}
                    }}
                }}
                document.getElementById("amountTable").innerText = count + " found";
            }}
        </script>
    <body class="container mt-4">
        <h2 class="text-center">DNA Viewer</h2>
        <div class="row">
                <class ="card card-body mb-4" >
                    <div class="row mb-4">
                        {plot_html}
                    </div>
                </div>
            <div class="col-md-12 mt-4">
            
            
            
            <div class="row d-flex align-items-center g-3">
                    <div class="col">
                        <input type="text" id="searchInput" class="form-control" placeholder="Search..." onkeyup="filterTable()">
                    </div>
                    <div class="col-auto">
                        <button class="btn btn-warning" type="button" id="amountTable"></button>
                    </div>
                </div>
            
            
            
                {table_html}
            </div>
        </div>
    </body>
    </html>
    """
    with open(f"{args.output}.html", 'w', encoding="utf-8") as file:
        file.write(html_template)



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

