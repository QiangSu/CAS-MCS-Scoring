#!/usr/bin/env python3
"""
Integrated standalone Python script for single-cell RNA-seq analysis using Scanpy.
This script combines multiple analysis and annotation strategies into one comprehensive
command-line workflow.

It performs the following steps:
a. Loads 10x Genomics data.
b. Performs detailed quality control (QC) and filtering.
c. Normalizes, log-transforms, and finds highly variable genes.
d. Scales data, runs PCA, computes a neighborhood graph, and performs Leiden clustering.
e. Calculates a Silhouette Score to evaluate clustering quality.
f. Generates a UMAP embedding for visualization.
g. Annotates cell types using a CellTypist model.
h. Performs marker gene analysis for different annotation levels.
i. Applies manual and marker-database-driven annotation strategies.
j. Calculates and exports a Cluster Annotation Score (CAS) for cell-type stability.
k. Calculates and exports a Marker Concordance Score (MCS) for each cell type.
l. Exports all annotations, marker lists, plots, and the final AnnData object.
"""

# --- Step 0: Imports ---
import os
import sys
import random
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import anndata
import celltypist
from celltypist import models
from sklearn.metrics import silhouette_score
import argparse

# ==============================================================================
# === HELPER FUNCTIONS ===
# ==============================================================================

def set_seed(seed):
    """Set random seeds for reproducibility."""
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    print(f"[INFO] Global random seed set to {seed}")

# ==============================================================================
# === MAIN ANALYSIS PIPELINE ===
# ==============================================================================

def main(args):
    """Main function to run the entire analysis pipeline."""
    # --- 1. Setup ---
    print("--- Initializing Integrated Analysis ---")
    set_seed(args.seed)
    sc.settings.verbosity = 3
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False, dpi_save=args.fig_dpi)

    os.makedirs(args.output_dir, exist_ok=True)
    sc.settings.figdir = args.output_dir

    print(f"Scanpy version: {sc.__version__}")
    print(f"Output directory: {os.path.abspath(args.output_dir)}")

    # --- 2. Load Data ---
    print("\n--- Step 2: Loading Data ---")
    adata = sc.read_10x_mtx(args.data_dir, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # --- 3. QC & Filtering ---
    print("\n--- Step 3: Quality Control and Filtering ---")
    adata.var['mt'] = adata.var_names.str.startswith(args.mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show=False)
    plt.savefig(os.path.join(args.output_dir, f"{args.prefix}_qc_violin_before_filtering.png"))
    plt.close()

    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_cells(adata, max_genes=args.max_genes)
    adata = adata[adata.obs.pct_counts_mt < args.max_mt_pct, :]
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"Filtered dims: {adata.n_obs} cells, {adata.n_vars} genes")

    # --- 4. Normalization, HVG, and Scaling ---
    print("\n--- Step 4: Normalization, HVG, Scaling ---")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvgs, flavor='seurat_v3')
    sc.pl.highly_variable_genes(adata, save=f"_{args.prefix}_hvg_plot.png", show=False)
    plt.close()
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)

    # --- 5. PCA, Clustering, UMAP & QC Score ---
    print("\n--- Step 5: Dimensionality Reduction and Clustering ---")
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50, random_state=args.seed)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save=f"_{args.prefix}_pca_variance.png", show=False)
    plt.close()

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=args.n_pcs)
    sc.tl.leiden(adata, resolution=args.leiden_res, random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)

    silhouette_avg = silhouette_score(adata.obsm['X_pca'][:, :args.n_pcs], adata.obs['leiden'])
    print(f"Average Silhouette Score for Leiden clustering: {silhouette_avg:.3f}")
    sc.pl.umap(
        adata, color='leiden', legend_loc='on data',
        title=f'Leiden Clusters (res={args.leiden_res}, nPCs={args.n_pcs})\nSilhouette Score: {silhouette_avg:.3f}',
        save=f"_{args.prefix}_umap_leiden.png", show=False
    )
    plt.close()

    # --- 6. CellTypist Annotation ---
    print("\n--- Step 6: CellTypist Annotation ---")
    model = models.Model.load(args.celltypist_model)
    predictions = celltypist.annotate(adata, model=model, majority_voting=True)
    
    adata.obs['ctpt_per_cell'] = predictions.predicted_labels['predicted_labels'].astype('category')
    adata.obs['ctpt_majority_vote'] = predictions.predicted_labels['majority_voting'].astype('category')
    if 'conf_score' in predictions.predicted_labels.columns:
        adata.obs['ctpt_confidence'] = predictions.predicted_labels['conf_score']

    sc.pl.umap(adata, color='ctpt_majority_vote', legend_loc='right margin', title='CellTypist Majority Voting Annotation',
                save=f"_{args.prefix}_umap_celltypist.png", show=False)
    plt.close()

    # --- 7. Manual and Marker-based Annotation ---
    print("\n--- Step 7: Applying Other Annotation Strategies ---")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=True, key_added="wilcoxon_leiden")
    
    # Manual Annotation from file
    if args.manual_map:
        try:
            map_df = pd.read_csv(args.manual_map)
            map_dict = dict(zip(map_df.iloc[:, 0].astype(str), map_df.iloc[:, 1].astype(str)))
            adata.obs['anno_manual'] = adata.obs['leiden'].map(map_dict).astype('category').fillna('Unknown')
            print("Manual annotation loaded and applied successfully.")
        except Exception as e:
            print(f"[WARNING] Could not apply manual annotation map: {e}")
            adata.obs['anno_manual'] = 'not_provided'
    else:
        adata.obs['anno_manual'] = 'not_provided'

    # Marker Concordance Score (MCS) and Annotation
    if args.marker_db:
        try:
            db_df = pd.read_csv(args.marker_db)
            db_markers_dict = {
                row.iloc[0]: set(str(m).strip().upper() for m in str(row.iloc[1]).split(','))
                for _, row in db_df.iterrows()
            }
            
            marker_names = adata.uns['wilcoxon_leiden']['names']
            mcs_results = []
            auto_map = {}

            for cluster_id in adata.obs['leiden'].cat.categories:
                cluster_markers = set(str(g).upper() for g in marker_names[cluster_id][:args.n_markers_mcs])
                concordance_scores = {"leiden_cluster": cluster_id}
                best_match, best_score = 'Unknown', 0.0

                for cell_type, db_markers in db_markers_dict.items():
                    if not db_markers: continue
                    score = len(cluster_markers.intersection(db_markers)) / len(db_markers)
                    concordance_scores[cell_type] = score
                    if score > best_score:
                        best_score, best_match = score, cell_type
                
                mcs_results.append(concordance_scores)
                auto_map[cluster_id] = best_match if best_score >= args.min_mcs_score else 'Unknown'

            adata.obs['anno_marker_based'] = adata.obs['leiden'].map(auto_map).astype('category')
            mcs_df = pd.DataFrame(mcs_results).set_index("leiden_cluster")
            
            mcs_path = os.path.join(args.output_dir, f"{args.prefix}_marker_concordance_score.csv")
            mcs_df.to_csv(mcs_path)
            print(f"Marker Concordance Score (MCS) matrix saved to: {mcs_path}")
        except Exception as e:
            print(f"[WARNING] Could not perform marker-based annotation. Error: {e}")
            adata.obs['anno_marker_based'] = 'not_provided'
    else:
        adata.obs['anno_marker_based'] = 'not_provided'

    # --- 8. Cluster Annotation Score (CAS) ---
    print("\n--- Step 8: Calculating Cluster Annotation Score (CAS) ---")
    cas_results = []
    for cluster_label, group_df in adata.obs.groupby('ctpt_majority_vote'):
        total = len(group_df)
        matching = (group_df['ctpt_per_cell'] == cluster_label).sum()
        cas = (matching / total) * 100 if total > 0 else 0
        cas_results.append({
            "Consensus_Label": cluster_label,
            "Total_Cells": total,
            "Matching_Per_Cell_Labels": matching,
            "CAS_Percentage": cas
        })
    cas_df = pd.DataFrame(cas_results).sort_values(by="Total_Cells", ascending=False)
    
    cas_path = os.path.join(args.output_dir, f"{args.prefix}_cluster_annotation_score.csv")
    cas_df.to_csv(cas_path, index=False)
    print(f"Cluster Annotation Score (CAS) summary saved to: {cas_path}")

    # --- 9. Export Results ---
    print("\n--- Step 9: Exporting Final Results ---")
    # Marker gene lists
    sc.tl.rank_genes_groups(adata, 'ctpt_majority_vote', method='wilcoxon', use_raw=True, key_added="wilcoxon_final")
    markers_df = sc.get.rank_genes_groups_df(adata, group=None, key="wilcoxon_final")
    markers_path = os.path.join(args.output_dir, f"{args.prefix}_marker_genes.csv")
    markers_df.to_csv(markers_path, index=False)
    print(f"Final marker genes saved to: {markers_path}")

    # Final marker heatmap
    sc.pl.heatmap(adata, var_names=markers_df['names'].unique()[:100],
                  groupby='ctpt_majority_vote', cmap='viridis', dendrogram=True, use_raw=True, show=False,
                  save=f"_{args.prefix}_final_marker_heatmap.png")
    plt.close()

    # All annotations CSV
    cols_to_save = ['leiden', 'ctpt_per_cell', 'ctpt_majority_vote', 'anno_manual', 'anno_marker_based']
    if 'ctpt_confidence' in adata.obs.columns: cols_to_save.insert(2, 'ctpt_confidence')
    annotations_path = os.path.join(args.output_dir, f"{args.prefix}_all_annotations.csv")
    adata.obs[[c for c in cols_to_save if c in adata.obs.columns]].to_csv(annotations_path)
    print(f"All cell annotations saved to: {annotations_path}")

    # Final AnnData object
    final_adata_path = os.path.join(args.output_dir, f"{args.prefix}_final_processed.h5ad")
    adata.write(final_adata_path)
    print(f"Final AnnData object saved to: {final_adata_path}")

    print("\n--- Integrated Pipeline Finished Successfully! ---")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Integrated Scanpy analysis pipeline with CAS and MCS scoring.")

    io = parser.add_argument_group('Input/Output')
    io.add_argument('--data_dir', required=True, type=str, help='Path to 10x data directory.')
    io.add_argument('--output_dir', required=True, type=str, help='Directory to save all outputs.')
    io.add_argument('--prefix', default='scanpy_integrated', type=str, help='Prefix for all output files.')

    qc = parser.add_argument_group('Quality Control')
    qc.add_argument('--min_genes', default=200, type=int, help='Min genes per cell.')
    qc.add_argument('--max_genes', default=7000, type=int, help='Max genes per cell.')
    qc.add_argument('--max_mt_pct', default=10.0, type=float, help='Max mitochondrial %% per cell.')
    qc.add_argument('--min_cells', default=3, type=int, help='Min cells per gene.')
    qc.add_argument('--mito_prefix', default='mt-', type=str, help="Prefix for mitochondrial genes ('mt-' or 'MT-').")

    analysis = parser.add_argument_group('Analysis Parameters')
    analysis.add_argument('--n_hvgs', default=10000, type=int, help='Number of highly variable genes.')
    analysis.add_argument('--n_pcs', default=8, type=int, help='Number of principal components to use.')
    analysis.add_argument('--leiden_res', default=0.9, type=float, help='Resolution for Leiden clustering.')
    
    annotation = parser.add_argument_group('Annotation Parameters')
    annotation.add_argument('--celltypist_model', required=True, type=str, help='Path to CellTypist .pkl model.')
    annotation.add_argument('--manual_map', type=str, default=None, help='Optional path to a 2-column CSV (cluster,cell_type) for manual annotation.')
    annotation.add_argument('--marker_db', type=str, default=None, help='Optional path to a 2-column CSV (cell_type,markers) for MCS calculation.')
    annotation.add_argument('--n_markers_mcs', default=25, type=int, help='Top N marker genes per cluster to use for MCS calculation.')
    annotation.add_argument('--min_mcs_score', default=0.05, type=float, help='Minimum MCS score to assign a cell type label.')

    other = parser.add_argument_group('Other Settings')
    other.add_argument('--seed', default=42, type=int, help='Random seed for reproducibility.')
    other.add_argument('--fig_dpi', default=300, type=int, help='Resolution (DPI) for saved figures.')

    args = parser.parse_args()
    main(args)