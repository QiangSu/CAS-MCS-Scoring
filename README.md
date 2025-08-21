# Post-Quantification Single-Cell Analysis using Scanpy and CellTypist
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository provides a suite of command-line pipelines for the analysis and cell type annotation of single-cell RNA-seq data. The workflows are built using [Scanpy](https://scanpy.readthedocs.io/) and leverage [CellTypist](https://www.celltypist.org/) for automated cell type annotation.

The provided scripts cover single-sample analysis, multi-sample integration with batch correction (Harmony), differential gene expression (DGE), and post-analysis quality control.

## Table of Contents
- [Available Workflows](#available-workflows)
- [Setup Instructions](#setup-instructions)
- [Usage](#usage)
  - [1. Single-Sample Analysis (Majority Voting)](#1-single-sample-analysis-majority-voting)
  - [2. Single-Sample Analysis (Per-Cell)](#2-single-sample-analysis-per-cell)
  - [3. Multi-Sample Integration & DGE](#3-multi-sample-integration--dge)
  - [4. Post-Analysis: Cluster Purity Calculation](#4-post-analysis-cluster-purity-calculation)
- [Citation](#citation)
- [License](#license)

## Available Workflows

This repository contains three core analysis pipelines and one post-processing utility script.

### 1. Single-Sample Analysis Pipelines
For end-to-end analysis of a **single 10x Genomics dataset**. They share a common workflow (QC, normalization, clustering, UMAP) but differ in their annotation strategy.

-   **`scripts/scanpy_pipeline_majority_voting.py`**: Performs clustering-based annotation. All cells within a given Leiden cluster are assigned the same consensus cell type label via a "majority vote".
-   **`scripts/scanpy_pipeline_per_cell.py`**: Performs per-cell annotation. Each cell is assigned an independent cell type label, which is useful for exploring cellular heterogeneity *within* clusters.

### 2. Multi-Sample Integration & DGE Pipeline
An advanced workflow for analyzing **two or more datasets**, such as comparing control vs. treated samples.

-   **`scripts/run_integration_analysis.py`**: Integrates multiple samples using Harmony to correct for batch effects. It performs a combined analysis including clustering, annotation, and **differential gene expression (DGE)** tests between specified conditions.

### 3. Post-Analysis Utilities
A script designed to be run after a primary analysis pipeline is complete.

-   **`scripts/calculate_cluster_purity.py`**: Evaluates annotation consistency from the **per-cell pipeline**. It calculates a "purity score" by measuring the percentage of cells within a cluster whose individual label matches the final consensus label for that cluster.

## Setup Instructions

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/miniconda.html) package manager.

### 1. Clone the Repository
```bash
git clone https://github.com/QiangSu/scanpy-celltypist-pipelines.git
cd scanpy-celltypist-pipelines
```

## Usage

All pipelines are run from the command line from the project’s root directory.

```bash
python scripts/calculate_cluster_purity.py \
    --input_file results/per_cell_output/WT_sample_annotations_per_cell_raw.csv \
    --output_file results/per_cell_output/WT_sample_cluster_purity_summary.csv
```

This will generate a file with the following columns:
- `Consensus_Cluster_Label`
- `Total_Cells_in_Cluster`
- `Cells_with_Matching_Individual_Label`
- `Purity_Percentage`


## 4. Two-Sample Integration and DGE Pipeline

This repository includes a comprehensive pipeline, `run_integration_analysis.py`, for the analysis of two or more scRNA-seq samples. The workflow integrates data using Harmony, performs clustering and several layers of annotation, and conducts differential gene expression (DGE) analysis between user-specified conditions.

### Configuration

This pipeline is controlled entirely via command-line arguments and two required configuration files:

**1. Sample Sheet (e.g., `sample_sheet.csv`)**
A CSV file that tells the script where to find each sample's data. It must contain `sample_id` and `path` columns.

*Example `sample_sheet.csv`:*
```csv
sample_id,path
WT,path/to/your/WT/filtered_matrices/
Treated,path/to/your/Treated/filtered_matrices/
```

**2. Manual Annotation Map (e.g., `manual_annotation_map.csv`)**
A CSV file for mapping Leiden cluster numbers to meaningful biological cell types. It must contain `leiden_cluster` and `cell_type_name` columns.

*Example `manual_annotation_map.csv`:*
```csv
leiden_cluster,cell_type_name
0,Astrocytes
1,Excitatory Neurons
2,Oligodendrocytes
```

### Usage Example

Once your configuration files are ready, you can run the pipeline from the project's root directory.

```bash
python scripts/run_integration_analysis.py \
    --sample_sheet sample_sheet.csv \
    --manual_annotation_map manual_annotation_map.csv \
    --celltypist_model models/Mouse_Whole_Brain.pkl \
    --output_dir results/integration_output \
    --output_prefix WT_vs_Treated_integrated \
    --dge_condition Treated \
    --dge_reference WT \
    --n_pcs 8 \
    --n_hvgs 10000
```
Use `python scripts/run_integration_analysis.py --help` to see all available options.


## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
