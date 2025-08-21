#!/usr/bin/env python3
"""
Standalone Python script for single-cell RNA-seq analysis using Scanpy.
This pipeline uses CellTypist with majority_voting=False for per-cell annotation.
Includes settings for reproducibility.

This script is designed to be run from the command line, with paths and parameters
provided as arguments.

It performs the following steps:
a. Loads 10x Genomics data.
b. Performs quality control (QC) and filtering.
c. Normalizes and log-transforms the data.
d. Identifies highly variable genes (HVGs).
e. Scales the data, runs PCA, computes neighbors, and performs Leiden clustering.
f. Generates a UMAP embedding for visualization.
g. Annotates cell types on a per-cell basis using a CellTypist model.
h. Exports annotations to CSV and saves the final AnnData object.
"""
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import celltypist
from celltypist import models

# <-- CHANGE: Standardized set_seed function for consistency
def set_seed(seed):
    """Set random seeds for reproducibility."""
    import random
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    print(f"Global random seed set to {seed}")

def main(args):
    """Main function to run the Scanpy analysis pipeline."""

    # --- 1. Reproducibility & Setup ---
    print("--- Initializing Analysis (Per-Cell Pipeline) ---")
    set_seed(args.seed) # Call the standardized seed function

    sc.settings.verbosity = 3
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False)

    os.makedirs(args.output_dir, exist_ok=True)
    sc.settings.figdir = args.output_dir

    print(f"Using Scanpy version {sc.__version__}")
    print(f"Data directory: {os.path.abspath(args.data_dir)}")
    print(f"Output directory: {os.path.abspath(args.output_dir)}")

    # --- 2. Load Data ---
    print("\n--- Step a: Loading Data ---")
    adata = sc.read_10x_mtx(args.data_dir, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    print(f"Loaded AnnData object: {adata.n_obs} cells x {adata.n_vars} genes")

    # --- 3. Quality Control and Filtering ---
    print("\n--- Step b: Quality Control and Filtering ---")
    adata.var['mt'] = adata.var_names.str.startswith(args.mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    n_cells_orig, n_genes_orig = adata.shape
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_cells(adata, max_genes=args.max_genes)
    adata = adata[adata.obs.pct_counts_mt < args.max_mt_pct, :]
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"Filtered dimensions: {adata.n_obs} cells ({adata.n_obs/n_cells_orig*100:.2f}%), {adata.n_vars} genes ({adata.n_vars/n_genes_orig*100:.2f}%)")

    # --- 4. Normalization & HVG Selection ---
    print("\n--- Step c/d: Normalization and HVG Selection ---")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvgs, flavor='seurat_v3')
    adata = adata[:, adata.var.highly_variable]
    print(f"Data subset to {adata.n_vars} highly variable genes.")

    # --- 5. Downstream Analysis (Scaling, PCA, Clustering, UMAP) ---
    print("\n--- Step e-h: Scaling, PCA, Clustering, and UMAP ---")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50, random_state=args.seed)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=args.n_pcs)
    sc.tl.leiden(adata, resolution=args.leiden_res, random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)
    print("Downstream analysis steps complete.")

    # --- 6. CellTypist Annotation (Per-Cell) ---
    print("\n--- Step i: CellTypist Annotation (majority_voting=False) ---")
    if not os.path.exists(args.model_path):
        raise FileNotFoundError(f"CellTypist model not found at: {args.model_path}")
    model = models.Model.load(args.model_path)
    predictions = celltypist.annotate(adata, model=model, majority_voting=False, mode='best match')
    
    adata.obs['ctpt_label'] = predictions.predicted_labels['predicted_labels'].astype('category')
    if 'conf_score' in predictions.predicted_labels.columns:
        adata.obs['ctpt_confidence'] = predictions.predicted_labels['conf_score']

    print("Per-cell annotation complete.")

    # --- 7. Export Results ---
    print("\n--- Step j: Exporting Results ---")

    # Save UMAP plot
    n_labels = len(adata.obs['ctpt_label'].cat.categories)
    sc.pl.umap(adata, color='ctpt_label', legend_loc='on data', 
              title=f'Per-Cell CellTypist Labels ({n_labels} types)',
              save=f"_{args.output_prefix}_umap_per_cell.png", show=False)
    plt.close()
    print(f"Saved UMAP plot to {args.output_dir}")

    # <-- CHANGE: Export a curated, final annotation table in addition to the raw one.
    print("Exporting curated per-cell annotations...")
    final_obs_cols = ['leiden', 'ctpt_label']
    if 'ctpt_confidence' in adata.obs.columns:
        final_obs_cols.append('ctpt_confidence')
    
    curated_csv_path = os.path.join(args.output_dir, f"{args.output_prefix}_per_cell_annotations.csv")
    adata.obs[final_obs_cols].to_csv(curated_csv_path)
    print(f"Curated per-cell annotations exported to: {curated_csv_path}")

    # Save final AnnData object with a standardized name
    adata_path = os.path.join(args.output_dir, f"{args.output_prefix}_final_annotated.h5ad")
    adata.write(adata_path)
    print(f"Final AnnData object saved to: {adata_path}")
    
    print("\n--- Per-Cell Pipeline finished successfully! ---")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Scanpy pipeline for scRNA-seq analysis with per-cell CellTypist annotation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    io_group = parser.add_argument_group('Input/Output')
    io_group.add_argument('--data_dir', type=str, required=True, 
                          help='Path to the directory containing 10x matrix, features, and barcodes files.')
    io_group.add_argument('--model_path', type=str, required=True, 
                          help='Path to the CellTypist .pkl model file.')
    io_group.add_argument('--output_dir', type=str, required=True, 
                          help='Path to the directory where all output files will be saved.')
    io_group.add_argument('--output_prefix', type=str, default='scanpy_analysis', 
                          help='Prefix for all output files.')
    
    qc_group = parser.add_argument_group('Quality Control Parameters')
    qc_group.add_argument('--mito_prefix', type=str, default='mt-', 
                          help="Prefix for mitochondrial genes (e.g., 'mt-' for mouse, 'MT-' for human).")
    qc_group.add_argument('--min_genes', type=int, default=200, 
                          help='Minimum number of genes detected per cell.')
    qc_group.add_argument('--max_genes', type=int, default=7000, 
                          help='Maximum number of genes detected per cell.')
    # <-- BUG FIX: Changed % to %% to prevent argparse crash on --help
    qc_group.add_argument('--max_mt_pct', type=float, default=10.0, 
                          help='Maximum percentage%% of mitochondrial counts allowed per cell.')
    qc_group.add_argument('--min_cells', type=int, default=3, 
                          help='Minimum number of cells a gene must be detected in.')

    analysis_group = parser.add_argument_group('Analysis Parameters')
    analysis_group.add_argument('--n_hvgs', type=int, default=10000, 
                                help='Number of highly variable genes to select for downstream analysis.')
    analysis_group.add_argument('--n_pcs', type=int, default=30, 
                                help='Number of principal components to use for neighborhood graph construction.')
    analysis_group.add_argument('--leiden_res', type=float, default=0.9, 
                                help='Resolution parameter for the Leiden clustering algorithm.')

    repro_group = parser.add_argument_group('Reproducibility')
    repro_group.add_argument('--seed', type=int, default=42, 
                             help='Random seed for PCA, Leiden, and UMAP to ensure deterministic results.')

    args = parser.parse_args()
    main(args)

