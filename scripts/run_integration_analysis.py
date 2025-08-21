#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Combined Scanpy analysis script for processing, integrating, and analyzing
two or more single-cell RNA-seq samples.

This script is designed as a portable command-line tool.

Workflow:
1.  Load datasets from a user-provided sample sheet.
2.  Concatenate and perform QC and filtering.
3.  Normalize, find HVGs, and scale.
4.  Perform PCA and integrate datasets using Harmony.
5.  Cluster with Leiden and generate UMAP.
6.  Perform automated cell type annotation with CellTypist.
7.  Find marker genes and apply a manual annotation map.
8.  Perform differential gene expression (DGE) between conditions.
9.  Generate key visualizations and export data.
"""

import os
import sys
import random
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import anndata
import harmonypy as hm
import celltypist
from celltypist import annotate, models

# --- Utility Functions (unchanged) ---
def extract_fraction_data(adata, output_dir, output_prefix, groupby_key, rank_genes_key, n_top_genes=5):
    """Calculates and saves dot plot fraction data."""
    # (The code for this function is exactly the same as in your original script)
    # ...
    pass # Placeholder for brevity, but include the full function here.

# --- Main Analysis Function ---
def main(args):
    """Main function to run the entire analysis pipeline."""
    
    # --- 1. Setup and Reproducibility ---
    print("--- Initializing Integration and DGE Pipeline ---")
    random.seed(args.seed)
    np.random.seed(args.seed)
    sc.settings.njobs = 1
    print(f"[INFO] Random seed set to: {args.seed}")

    sc.settings.verbosity = 3
    sc.logging.print_header()
    os.makedirs(args.output_dir, exist_ok=True)
    sc.settings.figdir = args.output_dir
    sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False, dpi_save=args.figure_dpi)
    print(f"Output Directory: {os.path.abspath(args.output_dir)}")
    print(f"Using {args.n_hvgs} HVGs and {args.n_pcs} PCs for integration.")

    # --- 2. Loading Datasets from Sample Sheet ---
    print("\n--- Step 1-2: Loading and Concatenating Datasets ---")
    try:
        sample_sheet = pd.read_csv(args.sample_sheet)
        assert 'sample_id' in sample_sheet.columns and 'path' in sample_sheet.columns
    except Exception as e:
        print(f"Error reading sample sheet '{args.sample_sheet}': {e}")
        sys.exit(1)

    adatas = {}
    for _, row in sample_sheet.iterrows():
        sample_id, path = row['sample_id'], row['path']
        print(f"... Loading sample: {sample_id} from {path}")
        adata = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)
        adata.var_names_make_unique()
        adata.obs['sample'] = sample_id
        adatas[sample_id] = adata
    
    adata = anndata.AnnData.concatenate(*adatas.values(), batch_key='sample', batch_categories=adatas.keys())
    print(f"Combined AnnData: {adata.n_obs} cells x {adata.n_vars} genes")

    # --- 3. QC and Filtering (using args) ---
    print("\n--- Step 3: Quality Control and Filtering ---")
    adata.var['mt'] = adata.var_names.str.startswith(args.mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_cells(adata, max_genes=args.max_genes)
    adata = adata[adata.obs.pct_counts_mt < args.max_mt_pct, :]
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"Filtered dimensions: {adata.n_obs} cells x {adata.n_vars} genes")

    # --- 4. Normalization, HVGs, Scaling ---
    print("\n--- Step 4: Normalization, HVG selection, Scaling ---")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvgs, batch_key='sample', flavor='seurat_v3')
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)

    # --- 5. PCA and Harmony Integration ---
    print("\n--- Step 5: PCA and Harmony Integration ---")
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50, random_state=args.seed)
    sc.external.pp.harmony_integrate(adata, key='sample', basis='X_pca', adjusted_basis='X_pca_harmony')
    
    # --- 6. Clustering and UMAP ---
    print("\n--- Step 6: Clustering and UMAP ---")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=args.n_pcs, use_rep='X_pca_harmony')
    sc.tl.leiden(adata, resolution=args.leiden_res, random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)
    sc.pl.umap(adata, color=['sample', 'leiden'], save=f"_{args.output_prefix}_umap_sample_leiden.png", show=False)
    plt.close()

    # --- 7. CellTypist Annotation ---
    print("\n--- Step 7: CellTypist Annotation ---")
    model = models.Model.load(args.celltypist_model)
    predictions = annotate(adata, model=model, majority_voting=True)
    adata.obs['cluster_ctpt'] = predictions.predicted_labels['majority_voting'].astype('category')
    sc.pl.umap(adata, color='cluster_ctpt', legend_loc='on data', save=f"_{args.output_prefix}_umap_celltypist.png", show=False)
    plt.close()

    # --- 8. Manual Annotation ---
    print("\n--- Step 8: Applying Manual Annotation Map ---")
    try:
        map_df = pd.read_csv(args.manual_annotation_map)
        # Ensure cluster column is string to match AnnData's leiden categories
        map_df['leiden_cluster'] = map_df['leiden_cluster'].astype(str)
        annotation_dict = pd.Series(map_df.cell_type_name.values, index=map_df.leiden_cluster).to_dict()
        adata.obs['cell_type_manual'] = adata.obs['leiden'].map(annotation_dict).astype('category')
        print("Manual annotation applied successfully.")
        sc.pl.umap(adata, color='cell_type_manual', legend_loc='on data', save=f"_{args.output_prefix}_umap_manual.png", show=False)
        plt.close()
    except Exception as e:
        print(f"Warning: Could not apply manual annotation map. Error: {e}")

    # --- 9. Differential Gene Expression (DGE) ---
    print("\n--- Step 9: Differential Gene Expression (DGE) Analysis ---")
    # Use the annotation specified by the user for DGE
    final_annotation_col = args.dge_annotation_column
    if final_annotation_col not in adata.obs.columns:
        print(f"Error: DGE annotation column '{final_annotation_col}' not found. Available: {list(adata.obs.columns)}")
        sys.exit(1)
    
    dge_results = {}
    for cell_type in adata.obs[final_annotation_col].cat.categories:
        print(f"... DGE for: {cell_type}")
        adata_celltype = adata[adata.obs[final_annotation_col] == cell_type].copy()
        if len(adata_celltype.obs['sample'].value_counts()) < 2:
            print(f"Skipping '{cell_type}': only one condition present.")
            continue

        sc.tl.rank_genes_groups(adata_celltype, groupby='sample', groups=[args.dge_condition], 
                                reference=args.dge_reference, method='wilcoxon', use_raw=True)
        result_df = sc.get.rank_genes_groups_df(adata_celltype, group=args.dge_condition)
        dge_results[cell_type] = result_df

    if dge_results:
        dge_path = os.path.join(args.output_dir, f"{args.output_prefix}_dge_results.xlsx")
        with pd.ExcelWriter(dge_path) as writer:
            for cell_type, df in dge_results.items():
                safe_name = ''.join(e for e in cell_type if e.isalnum())[:31]
                df.to_excel(writer, sheet_name=safe_name, index=False)
        print(f"DGE results saved to: {dge_path}")

    # --- 10. Save Final Object ---
    print("\n--- Step 10: Saving Final AnnData Object ---")
    final_adata_path = os.path.join(args.output_dir, f"{args.output_prefix}_final_processed.h5ad")
    adata.write(final_adata_path)
    print(f"Final AnnData object saved to: {final_adata_path}")
    print("\n--- ANALYSIS PIPELINE COMPLETE ---")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Scanpy integration and DGE analysis pipeline.")
    
    # --- I/O Arguments ---
    io_group = parser.add_argument_group('Input/Output')
    io_group.add_argument('--sample_sheet', type=str, required=True, help='Path to a CSV file with "sample_id" and "path" columns.')
    io_group.add_argument('--output_dir', type=str, required=True, help='Directory to save all outputs.')
    io_grop.add_argument('--output_prefix', type=str, default='scanpy_integration', help='Prefix for all output files.')
    io_group.add_argument('--celltypist_model', type=str, required=True, help='Path to CellTypist .pkl model.')
    io_group.add_argument('--manual_annotation_map', type=str, required=True, help='Path to a CSV with "leiden_cluster" and "cell_type_name" columns.')

    # --- QC Arguments ---
    qc_group = parser.add_argument_group('Quality Control')
    qc_group.add_argument('--min_genes', type=int, default=200, help='Min genes per cell.')
    qc_group.add_argument('--max_genes', type=int, default=7000, help='Max genes per cell.')
    qc_group.add_argument('--max_mt_pct', type=float, default=10.0, help='Max mitochondrial % per cell.')
    qc_group.add_argument('--min_cells', type=int, default=3, help='Min cells per gene.')
    qc_group.add_argument('--mito_prefix', type=str, default='mt-', help="Prefix for mitochondrial genes.")

    # --- Analysis Arguments ---
    an_group = parser.add_argument_group('Analysis Parameters')
    an_group.add_argument('--n_hvgs', type=int, default=10000, help='Number of highly variable genes.')
    an_group.add_argument('--n_pcs', type=int, default=8, help='Number of PCs to use for integration and clustering.')
    an_group.add_argument('--leiden_res', type=float, default=0.9, help='Resolution for Leiden clustering.')
    
    # --- DGE Arguments ---
    dge_group = parser.add_argument_group('Differential Gene Expression')
    dge_group.add_argument('--dge_annotation_column', type=str, default='cluster_ctpt', help='Column in adata.obs to use for grouping cells for DGE (e.g., cluster_ctpt, cell_type_manual).')
    dge_group.add_argument('--dge_condition', type=str, required=True, help='The condition/sample to test (e.g., "Treated"). Must match a name in sample_sheet.')
    dge_group.add_argument('--dge_reference', type=str, required=True, help='The reference condition (e.g., "WT"). Must match a name in sample_sheet.')

    # --- Other Arguments ---
    other_group = parser.add_argument_group('Other')
    other_group.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility.')
    other_group.add_argument('--figure_dpi', type=int, default=300, help='Resolution (DPI) for saved figures.')

    args = parser.parse_args()
    main(args)
