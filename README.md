# Post-Quantification Single-Cell Analysis using Scanpy and CellTypist

This repository contains two parallel Python pipelines for the automated analysis and cell type annotation of single-cell RNA-seq data, as described in [Your Paper/Preprint Title]. The analysis is performed using Scanpy, with cell type annotation handled by CellTypist.

## Pipelines

1.  **Majority Voting (`scanpy_pipeline_majority_voting.py`)**: Assigns a consensus cell type identity to all cells within a pre-computed Leiden cluster.
2.  **Per-Cell Annotation (`scanpy_pipeline_per_cell.py`)**: Assigns an independent cell type label to each individual cell.


## Setup

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/miniconda.html) package manager.

### 1. Clone the Repository
```bash
git clone https://github.com/your-username/my_scanpy_analysis_project.git
cd my_scanpy_analysis_project
```

### 2. Create the Conda Environment
This will install all the necessary packages with the correct versions.
```bash
conda env create -f environment/environment.yml
conda activate scanpy-env
```
*(The environment name, e.g., `scanpy-env`, will be defined inside your `environment.yml` file).*

### 3. Download Data and Models
The scripts require raw 10x Genomics count data and a pre-trained CellTypist model. Due to their size, these files are not included in the repository.

-   **Data**: Please follow the instructions in `data/README.md` to download and place the data.
-   **Model**: Please follow the instructions in `models/README.md` to download and place the model.

## Usage

After setting up the environment and downloading the required files, you can run the pipelines from the project root directory.

### Example: Running the Majority Voting Pipeline
```bash
python scripts/scanpy_pipeline_majority_voting.py \
    --data_dir data/filtered_feature_bc_matrix \
    --model_path models/Mouse_Whole_Brain.pkl \
    --output_dir results/majority_voting_output \
    --output_prefix WT_sample
```

### Example: Running the Per-Cell Pipeline
```bash
python scripts/scanpy_pipeline_per_cell.py \
    --data_dir data/filtered_feature_bc_matrix \
    --model_path models/Mouse_Whole_Brain.pkl \
    --output_dir results/per_cell_output \
    --output_prefix WT_sample
```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## 3. Post-Analysis: Cluster Purity Calculation

After running the per-cell annotation pipeline, you may want to quantify how consistent the individual cell labels are within each final cluster. The `calculate_cluster_purity.py` script is provided for this purpose.

It calculates a "purity score" for each consensus cluster label, defined as the percentage of cells within that cluster whose individual annotation matches the final consensus label. This script is most informative when run on the output of the **per-cell pipeline** (`scanpy_pipeline_per_cell.py`).

### Usage

The script takes the raw per-cell annotation CSV as input and produces a summary CSV as output.

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

