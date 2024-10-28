import json
from pathlib import Path
import argparse
import pandas as pd


def generate_demo_cell_marker_heatmap_options(cell_marker_all_path, species, tissue, selected_genes):
    # read xlsx here
    cell_marker_data = pd.read_excel(Path(__file__).parent / "Cell_marker_All.xlsx")
    selected_cell_marker_data = cell_marker_data[(cell_marker_data["species"] == species) & (cell_marker_data["tissue_class"] == tissue)]
    marker_genes_number = selected_cell_marker_data.groupby("cell_name")["marker"].count()
    targeted_rows = selected_cell_marker_data[selected_cell_marker_data["marker"].isin(selected_genes)]
    target_marker_genes_number = targeted_rows.groupby("cell_name")["marker"].count()

    target_marker_genes_proportion = (target_marker_genes_number / marker_genes_number).sort_values(ascending=False).dropna()
    x_headers = target_marker_genes_proportion.index.tolist()
    y_headers = targeted_rows.groupby("marker")["cell_name"].count().sort_values(ascending=False).index.tolist()

    # get data
    data = [[i, 0, f"{x:.2f}"] for i, x in enumerate(target_marker_genes_proportion)]
    for i in range(len(y_headers)):
        for j in range(len(x_headers)):
            data.append([j, i+1, "-" if targeted_rows[(targeted_rows["cell_name"] == x_headers[j]) & (targeted_rows["marker"] == y_headers[i])].empty else 1])
    return {"x_headers": x_headers, "y_headers": ["ALL"] + y_headers, "data": data}


if __name__ == '__main__':
    cell_marker_all_path = Path(__file__).parent / "Cell_marker_All.xlsx"

    parser = argparse.ArgumentParser()
    parser.add_argument('--gene_names_file', help='Selected gene names file path', required=True)
    parser.add_argument('--output_path', help='Output result directory path', required=True)
    args = parser.parse_args()

    gene_names_file = Path(args.gene_names_file)
    output_path = Path(args.output_path)

    with open(gene_names_file, 'r') as f:
        selected_genes = [line.strip() for line in f.readlines() if line.strip()]

    res_option = generate_demo_cell_marker_heatmap_options(cell_marker_all_path, "Mouse", "Brain", selected_genes)
    with open(output_path / "cell_marker_heatmap_options.json", 'w', encoding="UTF-8") as f:
        json.dump(res_option, f, ensure_ascii=False)
