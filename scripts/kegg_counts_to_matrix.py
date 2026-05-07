#!/usr/bin/env python3
"""
Convert KEGG counts H5 file to a filtered KEGG x genome matrix.

Reads the output from collate_globDB_kegg_annotations.py and filters based on:
- Minimum number of distinct KEGG annotations per genome
- Minimum number of distinct genomes per KEGG annotation

Outputs a KEGG ID x genome count matrix as H5 file.

Usage:
  python scripts/kegg_counts_to_matrix.py \
    --input output/globdb_kegg_counts.h5 \
    --output output/kegg_by_genome_matrix.h5 \
    --min-annotations 10 \
    --min-genomes 5
"""

from __future__ import annotations

import argparse
import os
from typing import Dict

import pandas as pd
import matplotlib.pyplot as plt


def read_kegg_counts(input_path: str) -> pd.DataFrame:
    """Read KEGG counts from HDF5 file."""
    return pd.read_hdf(input_path, key="kegg_counts")


def filter_by_thresholds(
    df: pd.DataFrame,
    min_annotations: int,
    min_genomes: int,
) -> pd.DataFrame:
    """Filter KEGG IDs and genomes based on minimum count thresholds.
    
    Args:
        df: DataFrame with columns [genome, kegg_id, count]
        min_annotations: Minimum distinct KEGG IDs per genome
        min_genomes: Minimum distinct genomes per KEGG ID
    
    Returns:
        Filtered DataFrame
    """
    # Filter genomes with at least min_annotations distinct KEGG IDs
    annotations_per_genome = df.groupby("genome")["kegg_id"].nunique()
    valid_genomes = annotations_per_genome[annotations_per_genome >= min_annotations].index
    df = df[df["genome"].isin(valid_genomes)]
    
    # Filter KEGG IDs found in at least min_genomes distinct genomes
    genomes_per_kegg = df.groupby("kegg_id")["genome"].nunique()
    valid_kegg = genomes_per_kegg[genomes_per_kegg >= min_genomes].index
    df = df[df["kegg_id"].isin(valid_kegg)]
    
    return df


def create_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Create KEGG x genome count matrix using pivot table.
    
    Args:
        df: Filtered DataFrame with columns [genome, kegg_id, count]
    
    Returns:
        Pivot table with KEGG IDs as rows, genomes as columns, counts as values
    """
    matrix = df.pivot_table(
        index="kegg_id",
        columns="genome",
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    return matrix


def write_matrix(matrix: pd.DataFrame, output_path: str) -> None:
    """Write count matrix to HDF5 file.
    
    Args:
        matrix: Count matrix (KEGG IDs x genomes)
        output_path: Path to output HDF5 file
    """
    print(f"Writing matrix ({matrix.shape[0]} KEGG IDs x {matrix.shape[1]} genomes) to {output_path}")
    os.makedirs(os.path.dirname(output_path), exist_ok=True) if os.path.dirname(output_path) else None
    matrix.to_hdf(output_path, key="kegg_matrix", mode="w")


def generate_histograms(df: pd.DataFrame, output_dir: str) -> None:
    """Generate and save histograms of annotation and genome distributions.
    
    Args:
        df: Filtered DataFrame with columns [genome, kegg_id, count]
        output_dir: Directory to save histogram images
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Histogram 1: Number of distinct KEGG annotations per genome
    annotations_per_genome = df.groupby("genome")["kegg_id"].nunique()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(annotations_per_genome, bins=50, edgecolor="black", alpha=0.7)
    ax.set_xlabel("Number of Distinct KEGG Annotations")
    ax.set_ylabel("Frequency (Number of Genomes)")
    ax.set_title("Distribution of KEGG Annotation Counts per Genome")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    hist1_path = os.path.join(output_dir, "annotations_per_genome_histogram.png")
    fig.savefig(hist1_path, dpi=300)
    print(f"Saved histogram to {hist1_path}")
    plt.close(fig)
    
    # Histogram 2: Number of distinct genomes per KEGG annotation
    genomes_per_kegg = df.groupby("kegg_id")["genome"].nunique()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(genomes_per_kegg, bins=50, edgecolor="black", alpha=0.7)
    ax.set_xlabel("Number of Distinct Genomes")
    ax.set_ylabel("Frequency (Number of KEGG Annotations)")
    ax.set_title("Distribution of Genome Counts per KEGG Annotation")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    hist2_path = os.path.join(output_dir, "genomes_per_kegg_histogram.png")
    fig.savefig(hist2_path, dpi=300)
    print(f"Saved histogram to {hist2_path}")
    plt.close(fig)


def build_arg_parser() -> argparse.ArgumentParser:
    """Build and return argument parser for command-line interface."""
    parser = argparse.ArgumentParser(
        description="Filter KEGG counts and create genome x KEGG annotation matrix."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input HDF5 file from collate_globDB_kegg_annotations.py"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output HDF5 file with filtered count matrix."
    )
    parser.add_argument(
        "-a", "--min-annotations",
        type=int,
        default=1,
        help="Minimum number of distinct KEGG annotations per genome (default: 1)."
    )
    parser.add_argument(
        "-g", "--min-genomes",
        type=int,
        default=1,
        help="Minimum number of distinct genomes per KEGG annotation (default: 1)."
    )
    return parser


def main() -> None:
    """Parse arguments, filter data, create matrix, and write output."""
    parser = build_arg_parser()
    args = parser.parse_args()

    print(f"Reading KEGG counts from {args.input}")
    df = read_kegg_counts(args.input)
    print(f"Loaded {len(df)} rows ({df['genome'].nunique()} genomes, {df['kegg_id'].nunique()} KEGG IDs)")

    print(f"Filtering: min_annotations={args.min_annotations}, min_genomes={args.min_genomes}")
    df_filtered = filter_by_thresholds(df, args.min_annotations, args.min_genomes)
    print(f"After filtering: {len(df_filtered)} rows ({df_filtered['genome'].nunique()} genomes, {df_filtered['kegg_id'].nunique()} KEGG IDs)")

    print("Creating matrix...")
    matrix = create_matrix(df_filtered)

    write_matrix(matrix, args.output)
    
    # Generate histograms from unfiltered data
    output_dir = os.path.dirname(args.output) or "output"
    generate_histograms(df, output_dir)


if __name__ == "__main__":
    main()
