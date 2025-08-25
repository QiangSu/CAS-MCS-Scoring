#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Comprehensive command-line pipeline for single-cell analysis and quantitative
annotation scoring (CAS & MCS).

This script performs an end-to-end analysis from 10x Genomics data, including
QC, clustering, and multiple annotation strategies. It calculates two key metrics:
- Cluster Annotation Score (CAS): The purity of CellTypist annotations within a
  consensus cluster.
- Marker Concordance Score (MCS): The average expression prevalence of a cluster's
  top 5 de novo marker genes, measuring internal consistency.
"""

import os
import sys
import argparse
import random
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import anndata
import celltypist
from celltypist import annotate, models
from sklearn.metrics import silhouette_score

# ==============================================================================
# === HELPER FUNCTIONS ===
# ==============================================================================

def reformat_dotplot_data(
    fraction_df: pd.DataFrame, top_genes_df: pd.DataFrame, output_dir: str,
    output_prefix: str, groupby_key: str
):
    """Reformats dot plot fraction data to a gene-centric sparse table."""
    print(f"[INFO] Reformatting dot plot data for '{groupby_key}'...")
    cell_types = top_genes_df['group'].unique().tolist()
    output_rows = []
    for _, row in top_genes_df.iterrows():
        gene, group = row['names'], row['group']
        fraction = fraction_df.loc[group, gene]
        new_row_data = {'Gene': gene, **{ct: '' for ct in cell_types}}
        new_row_data[group] = fraction
        output_rows.append(new_row_data)

    reformatted_df = pd.DataFrame(output_rows)[['Gene'] + cell_types]
    reformatted_csv_path = os.path.join(output_dir, f"{output_prefix}_dotplot_fractions_{groupby_key}_reformatted.csv")
    reformatted_df.to_csv(reformatted_csv_path, index=False)
    print(f"       -> Saved reformatted fraction data to: {reformatted_csv_path}")

def extract_fraction_data_and_calculate_mcs(
    adata: anndata.AnnData, output_dir: str, output_prefix: str,
    groupby_key: str, rank_genes_key: str, n_top_genes: int = 5
):
    """Calculates and saves expression fractions and the MCS."""
    print(f"[INFO] Calculating MCS and expression fractions for '{groupby_key}'...")
    if groupby_key not in adata.obs.columns or rank_genes_key not in adata.uns:
        print(f"[ERROR] Required keys not found for '{groupby_key}'. Skipping MCS calculation.")
        return

    marker_df = sc.get.rank_genes_groups_df(adata, key=rank_genes_key, group=None)
    top_genes_per_group = marker_df.groupby('group').head(n_top_genes)
    unique_top_genes = top_genes_per_group['names'].unique().tolist()
    
    # Calculate fraction of cells expressing genes
    data_df = sc.get.obs_df(adata, keys=[groupby_key] + unique_top_genes, use_raw=(adata.raw is not None))
    fraction_df = data_df.groupby(groupby_key).apply(lambda x: (x[unique_top_genes] > 0).mean())

    # --- MCS Calculation ---
    mcs_scores = {}
    for cell_type in top_genes_per_group['group'].unique():
        markers_for_type = top_genes_per_group[top_genes_per_group['group'] == cell_type]['names']
        # Extract the prevalence values for these specific markers for this cell type
        prevalence_values = fraction_df.loc[cell_type, markers_for_type]
        mcs_scores[cell_type] = prevalence_values.mean()

    mcs_df = pd.DataFrame.from_dict(mcs_scores, orient='index', columns=['MCS'])
    mcs_df.index.name = 'Cell_Type'
    mcs_csv_path = os.path.join(output_dir, f"{output_prefix}_marker_concordance_scores.csv")
    mcs_df.to_csv(mcs_csv_path)
    print(f"       -> Saved MCS scores to: {mcs_csv_path}")

    # --- Save full fraction data and reformat ---
    output_csv_path = os.path.join(output_dir, f"{output_prefix}_dotplot_fractions_{groupby_key}.csv")
    fraction_df.to_csv(output_csv_path)
    print(f"       -> Saved full fraction data to: {output_csv_path}")
    reformat_dotplot_data(fraction_df, top_genes_per_group, output_dir, output_prefix, groupby_key)

# ==============================================================================
# === MAIN ANALYSIS PIPELINE ===
# ==============================================================================

def main(args):
    """Main function to run the entire analysis pipeline."""
    print("--- Initializing CAS-MCS Scoring Pipeline ---")
    
    # --- Step 0: Setup and Reproducibility ---
    random.seed(args.seed)
    np.random.seed(args.seed)
    sc.settings.njobs = 1  # Enforce single-threaded behavior for determinism
    print(f"[INFO] Global random seed set to: {args.seed}")
    print("[INFO] Single-threaded execution enforced for reproducibility.")

    sc.settings.verbosity = 3
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False, dpi_save=args.fig_dpi)

    os.makedirs(args.output_dir, exist_ok=True)
    sc.settings.figdir = args.output_dir
    print(f"[INFO] Scanpy version: {sc.__version__}")
    print(f"[INFO] Output directory: {os.path.abspath(args.output_dir)}")

    # --- Step 1: Load Data ---
    print("\n--- Step 1: Loading Data ---")
    adata = sc.read_10x_mtx(args.data_dir, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    print(f"       -> Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # --- Step 2: QC & Filtering ---
    print("\n--- Step 2: Quality Control and Filtering ---")
    adata.var['mt'] = adata.var_names.str.startswith(args.mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Save pre-filtering QC plots
    fig_qc_before, axs = plt.subplots(1, 4, figsize=(18, 4))
    sc.pl.violin(adata, 'n_genes_by_counts', jitter=0.4, ax=axs[0], show=False)
    sc.pl.violin(adata, 'total_counts', jitter=0.4, ax=axs[1], show=False)
    sc.pl.violin(adata, 'pct_counts_mt', jitter=0.4, ax=axs[2], show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', ax=axs[3], show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, f"{args.prefix}_qc_plots_before_filtering.png"))
    plt.close()

    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_cells(adata, max_genes=args.max_genes)
    adata = adata[adata.obs.pct_counts_mt < args.max_pct_mt, :]
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"       -> Filtered dims: {adata.n_obs} cells, {adata.n_vars} genes")

    # --- Step 3: Normalization, HVG, and Scaling ---
    print("\n--- Step 3: Normalization, HVG, Scaling ---")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvgs, flavor='seurat_v3')
    sc.pl.highly_variable_genes(adata, save=f"_{args.prefix}_hvg_plot.png", show=False)
    plt.close()
    adata = adata[:, adata.var.highly_variable]

    sc.pp.scale(adata, max_value=10)

    # --- Step 4: PCA, Clustering, UMAP & QC Score ---
    print("\n--- Step 4: Dimensionality Reduction and Clustering ---")
    sc.tl.pca(adata, svd_solver='randomized', n_comps=50, random_state=args.seed)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save=f"_{args.prefix}_pca_variance.png", show=False)
    plt.close()

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=args.n_pcs)
    sc.tl.leiden(adata, resolution=args.resolution, random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)

    silhouette_avg = silhouette_score(adata.obsm['X_pca'][:, :args.n_pcs], adata.obs['leiden'])
    print(f"       -> Average Silhouette Score for Leiden clustering: {silhouette_avg:.3f}")
    sc.pl.umap(
        adata, color='leiden', legend_loc='on data',
        title=f'Leiden Clusters (res={args.resolution}, nPCs={args.n_pcs})\nSilhouette Score: {silhouette_avg:.3f}',
        save=f"_{args.prefix}_umap_leiden.png", show=False
    )
    plt.close()

    # --- Step 5: CellTypist Annotation (CAS) ---
    print("\n--- Step 5: CellTypist Annotation and CAS Calculation ---")
    model = models.Model.load(args.celltypist_model)
    predictions = celltypist.annotate(adata, model=model, majority_voting=not args.per_cell_annotation)
    
    label_col = 'predicted_labels' if args.per_cell_annotation else 'majority_voting'
    adata.obs['ctpt_individual_prediction'] = predictions.predicted_labels['predicted_labels'].astype('category')
    if 'conf_score' in predictions.predicted_labels.columns:
        adata.obs['ctpt_confidence'] = predictions.predicted_labels['conf_score']

    cluster2label = adata.obs.groupby('leiden')['ctpt_individual_prediction'].agg(lambda x: x.value_counts().idxmax()).to_dict()
    adata.obs['ctpt_consensus_prediction'] = adata.obs['leiden'].map(cluster2label).astype('category')
    
    sc.pl.umap(adata, color='ctpt_consensus_prediction', legend_loc='right margin', title='Cluster-Consensus CellTypist Annotation',
                save=f"_{args.prefix}_cluster_celltypist_umap.png", show=False)
    plt.close()

    # Calculate CAS (Purity)
    purity_results = []
    for cluster_name, group_df in adata.obs.groupby('ctpt_consensus_prediction'):
        total = len(group_df)
        matching = (group_df['ctpt_individual_prediction'] == cluster_name).sum()
        purity = (matching / total) * 100 if total > 0 else 0
        purity_results.append({
            "Consensus_Cell_Type": cluster_name,
            "Total_Cells_in_Cluster": total,
            "Matching_Individual_Predictions": matching,
            "Cluster_Annotation_Score_CAS (%)": purity
        })
    cas_df = pd.DataFrame(purity_results).sort_values(by="Total_Cells_in_Cluster", ascending=False)
    cas_output_path = os.path.join(args.output_dir, f"{args.prefix}_cluster_annotation_scores.csv")
    cas_df.to_csv(cas_output_path, index=False)
    print(f"       -> Saved CAS (Purity) scores to: {cas_output_path}")

    # --- Step 6: Marker Gene Analysis (MCS) ---
    print("\n--- Step 6: Marker Gene Analysis and MCS Calculation ---")
    marker_groupby_key = 'ctpt_consensus_prediction'
    marker_key = f"wilcoxon_{marker_groupby_key}"
    sc.tl.rank_genes_groups(adata, marker_groupby_key, method='wilcoxon', use_raw=True, key_added=marker_key)
    sc.pl.rank_genes_groups_dotplot(
        adata, n_genes=5, key=marker_key, groupby=marker_groupby_key,
        save=f"_{args.prefix}_markers_celltypist_dotplot.png", show=False
    )
    plt.close()
    
    extract_fraction_data_and_calculate_mcs(
        adata=adata, output_dir=args.output_dir, output_prefix=args.prefix,
        groupby_key=marker_groupby_key, rank_genes_key=marker_key, n_top_genes=5
    )

    # --- Step 7: Manual and Ratio-Based Annotation (Optional) ---
    print("\n--- Step 7: Applying Other Annotation Strategies ---")
    # Manual Annotation
    if args.manual_map_csv:
        try:
            map_df = pd.read_csv(args.manual_map_csv)
            manual_map = dict(zip(map_df['leiden_cluster'].astype(str), map_df['cell_type']))
            current_clusters = set(adata.obs['leiden'].cat.categories)
            missing_clusters = current_clusters - set(manual_map.keys())
            if missing_clusters:
                for cluster in missing_clusters: manual_map[cluster] = f'Unknown_{cluster}'
            adata.obs['cell_type_manual'] = adata.obs['leiden'].map(manual_map).astype('category')
            print("       -> Manual annotation counts:\n", adata.obs['cell_type_manual'].value_counts())
        except Exception as e:
            print(f"[WARNING] Could not apply manual annotation from {args.manual_map_csv}. Error: {e}")

    # Ratio-based Annotation
    if args.marker_db_csv:
        try:
            sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=True, key_added="wilcoxon_leiden")
            db_df = pd.read_csv(args.marker_db_csv)
            db_markers_dict = {row['Cell Type']: set(m.strip().upper() for m in row['Cell Marker'].split(',')) for _, row in db_df.iterrows() if isinstance(row['Cell Marker'], str)}
            marker_names = adata.uns['wilcoxon_leiden']['names']
            auto_map = {}
            for cluster_id in adata.obs['leiden'].cat.categories:
                cluster_markers = set(str(g).upper() for g in marker_names[cluster_id][:25])
                best_match, best_ratio = 'Unknown', 0.0
                for cell_type, db_markers in db_markers_dict.items():
                    if not db_markers: continue
                    ratio = len(cluster_markers.intersection(db_markers)) / len(db_markers)
                    if ratio > best_ratio: best_ratio, best_match = ratio, cell_type
                auto_map[cluster_id] = best_match if best_ratio >= 0.05 else 'Unknown (Low Score)'
            adata.obs['cell_type_auto_ratio'] = adata.obs['leiden'].map(auto_map).astype('category')
            print("       -> Ratio-based annotation counts:\n", adata.obs['cell_type_auto_ratio'].value_counts())
        except Exception as e:
            print(f"[WARNING] Could not perform ratio-based annotation. Error: {e}")

    # --- Step 8: Export Results ---
    print("\n--- Step 8: Exporting All Results ---")
    # Save final annotations CSV
    cols_to_save = ['leiden', 'ctpt_individual_prediction', 'ctpt_consensus_prediction']
    if 'ctpt_confidence' in adata.obs.columns: cols_to_save.insert(2, 'ctpt_confidence')
    if 'cell_type_manual' in adata.obs.columns: cols_to_save.append('cell_type_manual')
    if 'cell_type_auto_ratio' in adata.obs.columns: cols_to_save.append('cell_type_auto_ratio')
    
    annotations_path = os.path.join(args.output_dir, f"{args.prefix}_all_annotations.csv")
    adata.obs[[c for c in cols_to_save if c in adata.obs.columns]].to_csv(annotations_path)
    print(f"       -> All cell annotations saved to: {annotations_path}")

    # Save final AnnData object
    final_adata_path = os.path.join(args.output_dir, f"{args.prefix}_final_processed.h5ad")
    # Remove large, unnecessary items from .uns to keep file size down
    for key in list(adata.uns.keys()):
        if 'neighbors' in key or 'leiden' in key or 'pca' in key:
            del adata.uns[key]
    adata.write(final_adata_path)
    print(f"       -> Final AnnData object saved to: {final_adata_path}")

    print("\n--- CAS-MCS Scoring Pipeline Finished Successfully! ---")

# ==============================================================================
# === COMMAND-LINE INTERFACE SETUP ===
# ==============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
        Comprehensive command-line pipeline for single-cell analysis and quantitative
        annotation scoring (CAS & MCS).
        """
    )

    # --- Input/Output Arguments ---
    io_group = parser.add_argument_group('Input/Output')
    io_group.add_argument('--data_dir', type=str, required=True, help='Path to the directory...')
    io_group.add_argument('--output_dir', type=str, required=True, help='Path to the directory...')
    io_group.add_argument('--celltypist_model', type=str, required=True, help='Path to the model...')
    io_group.add_argument('--prefix', type=str, default='sc_analysis', help='Prefix for all output files.')
    io_group.add_argument('--marker_db_csv', type=str, default=None, help='(Optional) Path to a CSV file with canonical markers for ratio-based annotation.')
    io_group.add_argument('--manual_map_csv', type=str, default=None, help='(Optional) Path to a CSV file for manual cluster-to-celltype mapping. Must have "leiden_cluster" and "cell_type" columns.')
    
    # --- QC & Filtering Arguments ---
    qc_group = parser.add_argument_group('QC & Filtering')
    qc_group.add_argument('--min_genes', type=int, default=200, help='Minimum number of genes expressed in a cell for it to be kept.')
    qc_group.add_argument('--max_genes', type=int, default=7000, help='Maximum number of genes expressed in a cell for it to be kept.')
    qc_group.add_argument('--max_pct_mt', type=float, default=10.0, help='Maximum percentage of mitochondrial counts allowed in a cell.')
    qc_group.add_argument('--min_cells', type=int, default=3, help='Minimum number of cells a gene must be expressed in to be kept.')
    qc_group.add_argument('--mito_prefix', type=str, default='mt-', help='Prefix for mitochondrial genes (e.g., "mt-" for mouse, "MT-" for human).')

    # --- Analysis Arguments ---
    analysis_group = parser.add_argument_group('Analysis Parameters')
    analysis_group.add_argument('--n_hvgs', type=int, default=10000, help='Number of highly variable genes to select.')
    analysis_group.add_argument('--n_pcs', type=int, default=8, help='Number of principal components to use for neighborhood graph construction.')
    analysis_group.add_argument('--resolution', type=float, default=0.9, help='Resolution for Leiden clustering.')
    analysis_group.add_argument('--per_cell_annotation', action='store_true', help='If set, CellTypist performs per-cell annotation instead of majority voting. CAS will still be based on cluster consensus.')

    # --- Reproducibility & Other Arguments ---
    other_group = parser.add_argument_group('Other Settings')
    other_group.add_argument('--seed', default=42, type=int, help='Random seed for reproducibility.')
    other_group.add_argument('--fig_dpi', default=300, type=int, help='Resolution (DPI) for saved figures.')

    args = parser.parse_args()
    main(args)