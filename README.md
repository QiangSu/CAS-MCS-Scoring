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

### 1. Single-Sample Analysis (Majority Voting)

```bash
python scripts/scanpy_pipeline_majority_voting.py \
    --data_dir data/filtered_feature_bc_matrix \
    --model_path models/Mouse_Whole_Brain.pkl \
    --output_dir results/majority_voting_output \
    --output_prefix WT_sample
```

### 2. Single-Sample Analysis (Per-Cell)

```bash
python scripts/scanpy_pipeline_per_cell.py \
    --data_dir data/filtered_feature_bc_matrix \
    --model_path models/Mouse_Whole_Brain.pkl \
    --output_dir results/per_cell_output \
    --output_prefix WT_sample
```

### 3. Multi-Sample Integration & DGE
This pipeline requires two configuration files: a sample_sheet.csv and a manual_annotation_map.csv. See the example files in the root directory for the required format.

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
Use python scripts/run_integration_analysis.py --help to see all available options.

### 4. Post-Analysis: Cluster Purity Calculation
This script takes the raw annotation CSV from the per-cell pipeline as input.

```bash
python scripts/calculate_cluster_purity.py \
    --input_file results/per_cell_output/WT_sample_annotations_per_cell_raw.csv \
    --output_file results/per_cell_output/WT_sample_cluster_purity_summary.csv
```


---

## Data Simulation Pipeline

This repository also includes a two-step pipeline to generate realistic, ground-truth single-cell RNA-seq data. This is useful for benchmarking analysis methods or developing new tools.

### Step 1: Simulate Counts and Raw FASTQs (R)

The first script uses `splatter` to simulate a count matrix with known cell types and `polyester` to generate the corresponding raw paired-end FASTQ reads.

**Script:** `simulation/step1_simulate_counts_fastq.R`

**_Usage Example:_**
```bash
Rscript simulation/step1_simulate_counts_fastq.R \
    --ref_fasta path/to/your/Mus_musculus.GRCm38.cdna.all.fa \
    --output_dir results/simulation_step1_output \
    --n_cells 15000 \
    --n_groups 10

### Step 2: Add Barcodes and Create Cell Ranger FASTQs (Python)

The second script takes the raw FASTQs from Step 1, adds real 10x Genomics cell barcodes (from a whitelist) and random UMIs, and formats them into the R1/R2 file structure required by Cell Ranger.

**Script:** `simulation/step2_add_barcodes.py`

**_Usage Example:_**
```bash
python simulation/step2_add_barcodes.py \
    --input_dir results/simulation_step1_output/raw_fastq_files/ \
    --output_dir results/simulation_step2_cellranger_fastq/ \
    --whitelist path/to/your/3M-february-2018.txt.gz \
    --sample_name MySimulatedRun
```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
