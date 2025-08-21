#!/usr/bin/env python3
# ---
# Title: step2_add_barcodes.py
# Description: Takes simulated FASTQ files, adds real 10x barcodes and UMIs,
#              and structures them into Cell Ranger-compatible R1/R2 FASTQ files.
# Author: [Your Name]
# Date: [Current Date]
# ---
import argparse
import gzip
import random
import os
import sys
import glob

def generate_umi(length):
    """Generates a random UMI sequence."""
    return ''.join(random.choice('ATCG') for _ in range(length))

def main(args):
    print("--- Starting Step 2: Adding Barcodes/UMIs to create Cell Ranger FASTQs ---")
    os.makedirs(args.output_dir, exist_ok=True)

    # --- 1. Discover Input Files and Infer Cell Count ---
    input_files_pattern = os.path.join(args.input_dir, 'sample_*_1.fastq')
    input_files = sorted(glob.glob(input_files_pattern))
    n_cells = len(input_files)
    if n_cells == 0:
        print(f"ERROR: No input files found matching pattern: {input_files_pattern}")
        sys.exit(1)
    print(f"Found {n_cells} input cell FASTQ files to process.")

    # --- 2. Read Barcode Whitelist ---
    print(f"Reading barcode whitelist from: {args.whitelist}")
    valid_barcodes = []
    try:
        with gzip.open(args.whitelist, 'rt') as f:
            valid_barcodes = [line.strip() for line in f if line.strip()]
        print(f"Read {len(valid_barcodes)} barcodes from whitelist.")
        if len(valid_barcodes) < n_cells:
            raise ValueError(f"Whitelist has fewer barcodes ({len(valid_barcodes)}) than required cells ({n_cells}).")
    except Exception as e:
        print(f"ERROR: Failed to read or process barcode whitelist: {e}")
        sys.exit(1)

    # --- 3. Sample and Map Barcodes ---
    print(f"Sampling {n_cells} unique barcodes...")
    cell_barcodes_assigned = random.sample(valid_barcodes, n_cells)
    
    # --- 4. Process Files ---
    cr_r1_path = os.path.join(args.output_dir, f"{args.sample_name}_S1_L001_R1_001.fastq.gz")
    cr_r2_path = os.path.join(args.output_dir, f"{args.sample_name}_S1_L001_R2_001.fastq.gz")
    print(f"Writing combined output to:\n  R1: {cr_r1_path}\n  R2: {cr_r2_path}")

    total_reads_written = 0
    poly_t_tail_seq = 'T' * args.poly_t_len

    try:
        with gzip.open(cr_r1_path, 'wt') as r1_out, gzip.open(cr_r2_path, 'wt') as r2_out:
            for i, r1_path in enumerate(input_files):
                cell_barcode = cell_barcodes_assigned[i]
                
                if (i + 1) % 1000 == 0:
                    print(f"  Processing cell {i+1}/{n_cells}...")

                with open(r1_path, 'r') as r1_in:
                    while True:
                        header = r1_in.readline()
                        if not header: break
                        seq = r1_in.readline().strip()
                        plus = r1_in.readline()
                        qual = r1_in.readline().strip()

                        # Create Cell Ranger R1 (Barcode + UMI + PolyT)
                        umi = generate_umi(args.umi_len)
                        cr_r1_seq = cell_barcode + umi + poly_t_tail_seq
                        cr_r1_qual = 'I' * len(cr_r1_seq)

                        # Create Cell Ranger R2 (cDNA sequence)
                        cr_r2_seq = seq
                        cr_r2_qual = qual

                        # Write to combined gzipped files
                        r1_out.write(f"{header}{cr_r1_seq}\n+\n{cr_r1_qual}\n")
                        r2_out.write(f"{header}{cr_r2_seq}\n+\n{cr_r2_qual}\n")
                        total_reads_written += 1
    except Exception as e:
        print(f"FATAL ERROR during file processing: {e}")
        sys.exit(1)
        
    print("\n--- Processing Summary ---")
    print(f"Successfully processed {n_cells} input files.")
    print(f"Total reads written: {total_reads_written}")
    print("--- Simulation Step 2 Complete ---")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Add 10x barcodes/UMIs to simulated FASTQs.")
    parser.add_argument("--input_dir", type=str, required=True,
                        help="Directory containing raw FASTQ files from simulation step 1.")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory to save final Cell Ranger-compatible FASTQ files.")
    parser.add_argument("--whitelist", type=str, required=True,
                        help="Path to the gzipped 10x barcode whitelist file (e.g., 3M-february-2018.txt.gz).")
    parser.add_argument("--sample_name", type=str, default="SimulatedSample",
                        help="Sample name for the output Cell Ranger FASTQ files.")
    parser.add_argument("--barcode_len", type=int, default=16, help="Length of the cell barcode.")
    parser.add_argument("--umi_len", type=int, default=12, help="Length of the UMI.")
    parser.add_argument("--poly_t_len", type=int, default=10, help="Length of the poly-T tail to add to R1.")
    args = parser.parse_args()
    main(args)
