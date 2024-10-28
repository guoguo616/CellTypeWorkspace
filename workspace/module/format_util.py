import argparse
from pathlib import Path

import anndata
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite


def convert_h5ad_X_to_sparse(h5ad_file_path, output_dir):
    filename_no_ext = Path(h5ad_file_path).stem
    output_dir = Path(output_dir)
    adata = anndata.read_h5ad(h5ad_file_path)
    changed = False

    for key in adata.obsm.keys():
        if isinstance(adata.obsm[key], pd.DataFrame):
            print(f"Transforming: {key}: {type(adata.obsm[key])}")
            adata.obsm[key] = adata.obsm[key].values
            changed = True

    if changed:
        adata.write_h5ad(output_dir / f'{filename_no_ext}_sparse.h5ad')


def convert_h5ad_2_visium(h5ad_file_path, output_dir):
    adata = sc.read(h5ad_file_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    if not sparse.isspmatrix_csr(adata.X):
        adata.X = sparse.csr_matrix(adata.X)

    # mmwrite(output_dir / "expression_matrix.mtx", adata.X)
    adata.obs.to_csv(output_dir / 'cell_metadata.tsv', sep="\t")
    adata.var.to_csv(output_dir / 'gene_metadata.tsv', sep="\t")


if __name__ == '__main__':
    # my args should contain --method to tell the code which method to run
    # and other args to pass to the method
    parser = argparse.ArgumentParser()
    parser.add_argument('--method', help='Method to run', required=False)
    parser.add_argument('--h5ad_file_path', help='Input h5ad file path', required=False)
    parser.add_argument('--output_dir', help='Output directory path', required=False)
    args = parser.parse_args()

    # convert_h5ad_X_to_sparse(
    #     "D:\\st-network\\sample-data\\STANLY\\ExperimentalFolder\\h5ad\\sample-01.h5ad",
    #     "D:\\st-network\\cell-type\\cell-type-be\\workspace\\user_data\\local_test")

    output_path = Path("D:\\st-network\\cell-type\\cell-type-be\\workspace\\user_data\\local_test\\expression_matrix")
    convert_h5ad_2_visium(
        "D:\\st-network\\sample-data\\STANLY\\ExperimentalFolder\\h5ad\\sample-01.h5ad", output_path)

    if args.method == 'h5ad_X_2_sparse':
        convert_h5ad_X_to_sparse(args.h5ad_file_path, args.output_dir)
