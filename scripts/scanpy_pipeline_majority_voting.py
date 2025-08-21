#!/usr/bin/env python3
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
import argparse

def set_seed(seed):
    """Set random seeds for reproducibility."""
    import random
    import numpy as np
    import torch
    
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    # torch.manual_seed(seed) # Uncomment if using torch-dependent functions
    # torch.cuda.manual_seed_all(seed) # Uncomment if using GPU
    print(f"Global random seed set to {seed}")

def main(args):
    """Main function to run the Scanpy analysis pipeline."""
    # --- Setup ---
    print("--- Initializing Analysis (Majority Voting Pipeline) ---")
    set_seed(args.seed) # <-- CHANGE: Set seed for reproducibility
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
    
    # Filtering logic remains the same, using args for thresholds
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
    # <-- CHANGE: Added random_state for reproducibility
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50, random_state=args.seed)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=args.n_pcs)
    sc.tl.leiden(adata, resolution=args.leiden_res, random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)

    # --- Step i: CellTypist Annotation (Majority Voting) ---
    print("\n--- Step i: CellTypist Annotation (majority_voting=True) ---")
    if not os.path.exists(args.model_path):
        raise FileNotFoundError(f"CellTypist model not found at: {args.model_path}")
    model = models.Model.load(args.model_path)
    predictions = celltypist.annotate(adata, model=model, majority_voting=True, mode='best match')
    adata.obs['majority_voting'] = predictions.predicted_labels['majority_voting'].astype('category')
    
    # --- Step j: Export Results ---
    print("\n--- Step j: Exporting Results ---")
    # Save UMAP plot
    # <-- CHANGE: Standardized plot file naming
    sc.pl.umap(adata, color='majority_voting', legend_loc='on data', title='CellTypist Labels (Majority Voting)',
              save=f"_{args.output_prefix}_umap_majority_voting.png", show=False)
    plt.close()
    
    # <-- CHANGE: Export a clean, final consensus annotation table instead of the raw object
    print("Exporting final consensus annotations...")
    consensus_df = adata.obs[['leiden', 'majority_voting']].drop_duplicates().sort_values('leiden').set_index('leiden')
    csv_path = os.path.join(args.output_dir, f"{args.output_prefix}_consensus_annotations.csv")
    consensus_df.to_csv(csv_path)
    print(f"Consensus annotation table saved to: {csv_path}")
    
    # Save final object
    adata_path = os.path.join(args.output_dir, f"{args.output_prefix}_final_annotated.h5ad")
    adata.write(adata_path)
    print(f"Final AnnData object saved to: {adata_path}")
    print("\n--- Majority Voting Pipeline finished successfully! ---")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Scanpy pipeline with CellTypist (majority_voting=True).")
    
    # <-- CHANGE: Organized arguments into groups for cleaner --help output
    io_group = parser.add_argument_group('Input/Output')
    io_group.add_argument('--data_dir', type=str, required=True, help='Path to 10x data directory (containing matrix.mtx.gz, etc.).')
    io_group.add_argument('--model_path', type=str, required=True, help='Path to the CellTypist .pkl model file.')
    io_group.add_argument('--output_dir', type=str, required=True, help='Directory to save all outputs.')
    io_group.add_argument('--output_prefix', type=str, default='scanpy_run', help='Prefix for all output file names.')
    
    qc_group = parser.add_argument_group('Quality Control')
    qc_group.add_argument('--mito_prefix', type=str, default='mt-', help="Prefix for mitochondrial genes (e.g., 'mt-' or 'MT-').")
    qc_group.add_argument('--min_genes', type=int, default=200, help='Minimum number of genes expressed per cell.')
    qc_group.add_argument('--max_genes', type=int, default=7000, help='Maximum number of genes expressed per cell.')
    # <-- BUG FIX: Changed % to %% to prevent argparse error
    qc_group.add_argument('--max_mt_pct', type=float, default=10.0, help='Maximum allowed mitochondrial %% per cell.')
    qc_group.add_argument('--min_cells', type=int, default=3, help='Minimum number of cells a gene must be expressed in.')

    analysis_group = parser.add_argument_group('Analysis Parameters')
    analysis_group.add_argument('--n_hvgs', type=int, default=10000, help='Number of highly variable genes to select.')
    analysis_group.add_argument('--n_pcs', type=int, default=30, help='Number of principal components to use for neighborhood graph.')
    analysis_group.add_argument('--leiden_res', type=float, default=0.9, help='Resolution parameter for Leiden clustering.')

    repro_group = parser.add_argument_group('Reproducibility')
    # <-- CHANGE: Added seed argument
    repro_group.add_argument('--seed', type=int, default=42, help='Random seed for PCA, Leiden, and UMAP to ensure deterministic results.')

    args = parser.parse_args()
    main(args)
