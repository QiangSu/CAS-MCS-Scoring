#!/usr/bin/env python
"""
Standalone Python script for single-cell RNA-seq analysis using Scanpy.
This pipeline uses CellTypist with majority_voting=True.

It performs:
a. Data loading
b. QC and filtering
c. Normalization and HVG selection
d. Scaling, PCA, Neighborhood Graph, and Leiden Clustering
e. UMAP visualization
f. CellTypist annotation (cluster-level consensus)
g. Exporting results
"""
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os
import celltypist
from celltypist import models
from sklearn.metrics import silhouette_score
import argparse # <-- Import argparse

def main(args):
    """Main function to run the Scanpy analysis pipeline."""
    # --- Setup ---
    print("--- Initializing Analysis (Majority Voting Pipeline) ---")
    sc.settings.verbosity = 3
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False)

    os.makedirs(args.output_dir, exist_ok=True)
    sc.settings.figdir = args.output_dir

    print(f"Using Scanpy version {sc.__version__}")
    print(f"Data directory: {os.path.abspath(args.data_dir)}")
    print(f"Output directory: {os.path.abspath(args.output_dir)}")

    # --- Step a: Load Data ---
    print("\n--- Step a: Loading Data ---")
    adata = sc.read_10x_mtx(args.data_dir, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # --- Step b: QC and Filtering ---
    print("\n--- Step b: Quality Control and Filtering ---")
    adata.var['mt'] = adata.var_names.str.startswith(args.mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # (Filtering logic remains the same, using args for thresholds)
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_cells(adata, max_genes=args.max_genes)
    adata = adata[adata.obs.pct_counts_mt < args.max_mt_pct, :]
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"Filtered dimensions: {adata.n_obs} cells x {adata.n_vars} genes")

    # --- Step c & d: Normalization & HVG ---
    print("\n--- Step c/d: Normalization and HVG Selection ---")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvgs, flavor='seurat_v3')
    adata = adata[:, adata.var.highly_variable]

    # --- Step e-h: Scaling, PCA, Clustering, UMAP ---
    print("\n--- Step e-h: Downstream Analysis ---")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=args.n_pcs)
    sc.tl.leiden(adata, resolution=args.leiden_res)
    sc.tl.umap(adata)

    # --- Step i: CellTypist Annotation (Majority Voting) ---
    print("\n--- Step i: CellTypist Annotation (majority_voting=True) ---")
    if not os.path.exists(args.model_path):
        raise FileNotFoundError(f"CellTypist model not found at: {args.model_path}")
    model = models.Model.load(args.model_path)
    predictions = celltypist.annotate(adata, model=model, majority_voting=True, mode='best match')
    adata.obs['ctpt_label'] = predictions.predicted_labels['majority_voting'].astype('category')
    
    # --- Step j: Export Results ---
    print("\n--- Step j: Exporting Results ---")
    # Save UMAP plot
    sc.pl.umap(adata, color='ctpt_label', legend_loc='on data', title='CellTypist Labels (Majority Voting)',
              save=f"_{args.output_prefix}_celltypist_majority_voting.png", show=False)
    plt.close()
    
    # Export annotations
    csv_path = os.path.join(args.output_dir, f"{args.output_prefix}_annotations_majority_voting.csv")
    predictions.predicted_labels.to_csv(csv_path)
    
    # Save final object
    adata_path = os.path.join(args.output_dir, f"{args.output_prefix}_final_majority_voting.h5ad")
    adata.write(adata_path)
    print(f"Final AnnData object saved to: {adata_path}")
    print("\n--- Majority Voting Pipeline finished successfully! ---")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Scanpy pipeline with CellTypist (majority_voting=True).")
    
    # --- I/O Arguments ---
    parser.add_argument('--data_dir', type=str, required=True, help='Path to 10x data directory.')
    parser.add_argument('--model_path', type=str, required=True, help='Path to CellTypist .pkl model.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save outputs.')
    parser.add_argument('--output_prefix', type=str, default='scanpy_run', help='Prefix for output files.')
    
    # --- QC Arguments ---
    parser.add_argument('--mito_prefix', type=str, default='mt-', help="Prefix for mitochondrial genes ('mt-' or 'MT-').")
    parser.add_argument('--min_genes', type=int, default=200, help='Min genes per cell.')
    parser.add_argument('--max_genes', type=int, default=7000, help='Max genes per cell.')
    parser.add_argument('--max_mt_pct', type=float, default=10.0, help='Max mitochondrial % per cell.')
    parser.add_argument('--min_cells', type=int, default=3, help='Min cells per gene.')

    # --- Analysis Arguments ---
    parser.add_argument('--n_hvgs', type=int, default=10000, help='Number of highly variable genes.')
    parser.add_argument('--n_pcs', type=int, default=30, help='Number of principal components to use.')
    parser.add_argument('--leiden_res', type=float, default=0.9, help='Resolution for Leiden clustering.')

    args = parser.parse_args()
    main(args)
