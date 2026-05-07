#!/usr/bin/env python3
"""
Filter KEGG x genome matrix to a subset of KEGG genes.

Reads a count matrix and a list of target KEGG gene IDs, keeps only rows 
matching the target genes, and writes the filtered matrix.

Usage:
  python scripts/filter_kegg_matrix.py \
    --input output/kegg_by_genome_matrix.h5 \
    --genes data/kegg_genes.csv \
    --output output/kegg_subset_matrix.h5
"""

from __future__ import annotations

import argparse
import os
from typing import Set

import pandas as pd


def load_kegg_genes(genes_file: str) -> Set[str]:
    """Load KEGG gene IDs from CSV file (first column is the index).
    
    Args:
        genes_file: Path to CSV file with KEGG IDs in first column
    
    Returns:
        Set of KEGG gene IDs to filter
    """
    df = pd.read_csv(genes_file, index_col=0)
    # exclude NaN index values -- these occur when there are blanks in the CSV
    not_nan = df.index[~df.index.isna()] 
    gene_ids = not_nan.astype(str)
    print(f"Loaded {len(gene_ids)} KEGG genes from {genes_file}")
    return set(gene_ids)


def filter_matrix(matrix: pd.DataFrame, target_genes: Set[str]) -> pd.DataFrame:
    """Filter matrix rows to keep only those in target_genes.
    
    Args:
        matrix: KEGG x genome count matrix
        target_genes: Set of KEGG gene IDs to keep
    
    Returns:
        Filtered matrix with only target genes as rows
    """
    # Find intersection of matrix rows and target genes
    matrix_genes = set(matrix.index.astype(str))
    common_genes = matrix_genes & target_genes
    
    print(f"Matrix contains {len(matrix_genes)} genes")
    print(f"Found {len(common_genes)} matching genes ({100 * len(common_genes) / len(target_genes):.1f}% of target)")
    
    # Missing genes from target set
    missing = target_genes - matrix_genes
    if missing:
        print(f"Warning: {len(missing)} target genes not found in matrix")
        if len(missing) <= 20:
            print(f"  Missing: {', '.join(sorted(str(x) for x in missing))}")
    
    # Filter and return
    filtered = matrix.loc[matrix.index.astype(str).isin(common_genes)]
    return filtered


def write_filtered_matrix(matrix: pd.DataFrame, output_path: str) -> None:
    """Write filtered matrix to HDF5 file.
    
    Args:
        matrix: Filtered count matrix
        output_path: Path to output HDF5 file
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True) if os.path.dirname(output_path) else None
    print(f"Writing filtered matrix ({matrix.shape[0]} genes x {matrix.shape[1]} genomes) to {output_path}")
    matrix.to_hdf(output_path, key="kegg_matrix", mode="w")


def build_arg_parser() -> argparse.ArgumentParser:
    """Build and return argument parser for command-line interface."""
    parser = argparse.ArgumentParser(
        description="Filter KEGG x genome matrix to a subset of genes."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input HDF5 matrix file from kegg_counts_to_matrix.py"
    )
    parser.add_argument(
        "-g", "--genes",
        required=True,
        help="CSV file with KEGG gene IDs in first column"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output HDF5 file with filtered matrix"
    )
    return parser


def main() -> None:
    """Parse arguments, load genes, filter matrix, and write output."""
    parser = build_arg_parser()
    args = parser.parse_args()

    print(f"Reading matrix from {args.input}")
    matrix = pd.read_hdf(args.input, key="kegg_matrix")
    print(f"Loaded matrix: {matrix.shape[0]} genes x {matrix.shape[1]} genomes")
    
    target_genes = load_kegg_genes(args.genes)
    
    filtered_matrix = filter_matrix(matrix, target_genes)
    
    write_filtered_matrix(filtered_matrix, args.output)


if __name__ == "__main__":
    main()
