#!/usr/bin/env Rscript
# ---
# Title: step1_simulate_counts_fastq.R
# Description: Simulates scRNA-seq count data with Splatter and generates
#              corresponding FASTQ files with Polyester.
# Author: [Your Name]
# Date: [Current Date]
# ---

# --- Argument Parsing ---
suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description="Simulate scRNA-seq counts and FASTQ files.")
parser$add_argument("--ref_fasta", type="character", required=TRUE,
                    help="Path to the reference transcriptome FASTA file (e.g., GRCm38.cdna.all.fa).")
parser$add_argument("--output_dir", type="character", required=TRUE,
                    help="Base directory to save all simulation outputs.")
parser$add_argument("--n_cells", type="integer", default=15000,
                    help="Number of cells to simulate.")
parser$add_argument("--n_genes", type="integer", default=25000,
                    help="Number of genes (transcripts) to simulate.")
parser$add_argument("--n_groups", type="integer", default=10,
                    help="Number of distinct cell type groups to simulate.")
parser$add_argument("--read_len", type="integer", default=150,
                    help="Length of simulated reads.")
parser$add_argument("--seed", type="integer", default=123,
                    help="Random seed for reproducibility.")
args <- parser$parse_args()

# --- 0: Set Up ---
print("--- Starting Simulation Pipeline: Counts with Splatter and FASTQ with Polyester ---")
set.seed(args$seed)

# --- 1: Load Libraries ---
print("Loading required libraries...")
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(polyester))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(Matrix))

# --- 2: Define Paths from Arguments ---
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
ground_truth_matrix_file <- file.path(args$output_dir, "ground_truth_counts.csv")
ground_truth_celltypes_file <- file.path(args$output_dir, "ground_truth_celltypes.tsv")
output_fastq_dir <- file.path(args$output_dir, "raw_fastq_files")
dir.create(output_fastq_dir, showWarnings = FALSE)

# --- 3: Load Reference Transcript IDs ---
print(paste("Loading reference transcript IDs from:", args$ref_fasta))
if (!file.exists(args$ref_fasta)) {
    stop(paste("Reference FASTA file not found:", args$ref_fasta))
}
fasta_data <- readDNAStringSet(args$ref_fasta)
transcript_ids <- sapply(strsplit(names(fasta_data), " "), `[`, 1)
transcript_ids <- sub("\\..*", "", transcript_ids)
transcript_ids <- unique(transcript_ids)
print(paste("Loaded", length(transcript_ids), "unique transcript IDs."))

# --- 4: Simulate Counts with Splatter ---
print("Simulating single-cell RNA-seq counts with Splatter...")
n_genes_final <- min(args$n_genes, length(transcript_ids))
if (n_genes_final < args$n_genes) {
    warning(paste("Target n_genes", args$n_genes, "is > available unique transcripts.",
                  "Using", n_genes_final, "genes instead."))
}
selected_gene_ids <- sample(transcript_ids, n_genes_final, replace = FALSE)

params <- newSplatParams()
params <- setParam(params, "nGenes", n_genes_final)
params <- setParam(params, "batchCells", args$n_cells)
params <- setParam(params, "group.prob", rep(1/args$n_groups, args$n_groups))
params <- setParam(params, "de.prob", 0.8)
params <- setParam(params, "de.facLoc", 3)
params <- setParam(params, "de.facScale", 0.6)
params <- setParam(params, "bcv.common", 0.2)
params <- setParam(params, "dropout.type", "experiment")
params <- setParam(params, "seed", args$seed)
sim_data <- splatSimulate(params, method = "groups", verbose = TRUE)

counts <- counts(sim_data)
rownames(counts) <- selected_gene_ids
colnames(counts) <- paste0("Cell_", 1:args$n_cells)
sim_cell_types <- colData(sim_data)$Group
names(sim_cell_types) <- colnames(counts)

print(paste("Saving ground truth count matrix to:", ground_truth_matrix_file))
write.csv(as.matrix(counts), file = ground_truth_matrix_file, quote = FALSE)

print(paste("Saving ground truth cell types to:", ground_truth_celltypes_file))
write.table(data.frame(CellID=names(sim_cell_types), CellType=sim_cell_types),
            file = ground_truth_celltypes_file, sep="\t", row.names = FALSE, quote = FALSE)

# --- 5: Simulate FASTQ Files with Polyester ---
print("Simulating FASTQ files with Polyester...")
ref_transcriptome <- fasta_data
names(ref_transcriptome) <- sub("\\..*", "", sapply(strsplit(names(ref_transcriptome), " "), `[`, 1))
ref_transcriptome <- ref_transcriptome[!duplicated(names(ref_transcriptome))]

common_transcripts <- intersect(rownames(counts), names(ref_transcriptome))
print(paste("Found", length(common_transcripts), "common transcripts between counts and reference."))
counts_for_polyester <- counts[common_transcripts, ]
ref_transcriptome_filtered <- ref_transcriptome[common_transcripts]

temp_fasta_file <- file.path(args$output_dir, "temp_filtered_transcriptome.fa")
writeXStringSet(ref_transcriptome_filtered, filepath = temp_fasta_file, format = "fasta")

simulate_experiment_countmat(fasta = temp_fasta_file,
                             readmat = as.matrix(counts_for_polyester),
                             outdir = output_fastq_dir,
                             paired = TRUE,
                             readlen = args$read_len,
                             error_rate = 0.005,
                             seed = args$seed)

unlink(temp_fasta_file)
print(paste("FASTQ simulation complete. Raw FASTA files saved to:", output_fastq_dir))

# --- 6: Rename Polyester Output Files ---
print("Renaming Polyester output files (.fasta -> .fastq)...")
fasta_files <- list.files(output_fastq_dir, pattern = "\\.fasta$", full.names = TRUE)
if (length(fasta_files) > 0) {
    fastq_files <- sub("\\.fasta$", ".fastq", fasta_files)
    file.rename(from = fasta_files, to = fastq_files)
    print(paste("Successfully renamed", length(fasta_files), "files to .fastq format."))
} else {
    print("No .fasta files found to rename.")
}

print("--- Simulation Step 1 Complete ---")
