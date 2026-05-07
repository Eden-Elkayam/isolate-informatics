#!/usr/bin/env python3
"""
Evaluate pathway presence and completeness in genomes.

Takes a KEGG matrix (genes x genomes), pathway definitions with boolean
expressions for marker genes, and a complete gene list for completeness
calculation.

Presence evaluation: Boolean expressions in pathway_definitions.csv identify
  hallmark/marker genes. Uses AND/OR/NOT logic.

Completeness calculation: All genes per pathway from kegg_genes.csv.
  Calculates fraction of genes present (0-1).

Usage:
  python scripts/evaluate_pathway_logic.py \
    --matrix output/kegg_subset_matrix.h5 \
    --pathways data/pathway_definitions.csv \
    --genes data/kegg_genes.csv \
    --output output/pathway_presence_matrix.h5
"""

from __future__ import annotations

import argparse
import os
import re
from typing import Set

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm


def load_matrix(matrix_file: str) -> pd.DataFrame:
    """Load KEGG gene x genome matrix.
    
    Args:
        matrix_file: Path to HDF5 matrix file
    
    Returns:
        DataFrame with KEGG IDs as rows, genomes as columns
    """
    df = pd.read_hdf(matrix_file, key="kegg_matrix")
    print(f"Loaded matrix: {len(df)} genes x {len(df.columns)} genomes")
    return df


def load_pathway_definitions(pathways_file: str) -> pd.DataFrame:
    """Load pathway definitions with boolean expressions.
    
    Args:
        pathways_file: Path to CSV file with columns: pathway, expression, meaning
    
    Returns:
        DataFrame with pathway definitions
    """
    df = pd.read_csv(pathways_file)
    
    required_cols = ["pathway", "expression"]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in pathway file: {missing}")
    
    print(f"Loaded {len(df)} pathway definitions")
    return df


def load_kegg_genes_by_pathway(genes_file: str) -> dict[str, Set[str]]:
    """Load KEGG genes and group by pathway.
    
    Args:
        genes_file: Path to CSV file with KEGG genes and pathway column
    
    Returns:
        Dictionary mapping pathway name to set of KEGG IDs
    """
    df = pd.read_csv(genes_file, index_col=0)
    
    if "pathway" not in df.columns:
        raise ValueError(f"pathway column not found in {genes_file}")
    
    # Group genes by pathway
    pathway_genes = {}
    for pathway, group in df.groupby("pathway"):
        gene_ids = set(group.index.astype(str))
        pathway_genes[pathway] = gene_ids
    
    print(f"Loaded {len(df)} genes across {len(pathway_genes)} pathways")
    return pathway_genes


def extract_kegg_ids(expression: str) -> Set[str]:
    """Extract all KEGG IDs from a boolean expression.
    
    Args:
        expression: Boolean expression string
    
    Returns:
        Set of KEGG IDs found in the expression
    """
    # Find all K##### patterns
    kegg_pattern = r'K\d{5}'
    return set(re.findall(kegg_pattern, expression))


def create_gene_presence_dict(matrix: pd.DataFrame, genome: str) -> dict[str, bool]:
    """Create a dictionary of gene presence for a single genome.
    
    Args:
        matrix: KEGG matrix
        genome: Genome ID
    
    Returns:
        Dictionary mapping KEGG IDs to boolean presence (count > 0)
    """
    genome_data = matrix[genome]
    return {kegg_id: bool(count > 0) for kegg_id, count in genome_data.items()}


def evaluate_expression(expression: str, gene_presence: dict[str, bool]) -> bool:
    """Evaluate a boolean expression for a single genome.
    
    Args:
        expression: Boolean expression over KEGG IDs (e.g., "K01601 AND (K00855 OR K28462)")
        gene_presence: Dictionary mapping KEGG IDs to boolean presence
    
    Returns:
        Boolean result of evaluating the expression
    """
    # Replace KEGG IDs with their boolean values
    # Make a copy to work with
    eval_expr = expression
    
    # Find all KEGG IDs in the expression
    kegg_ids = extract_kegg_ids(expression)
    
    # Replace each KEGG ID with its boolean value
    for kegg_id in kegg_ids:
        # Get presence, defaulting to False if not in matrix
        presence = gene_presence.get(kegg_id, False)
        # Replace with Python boolean (case-sensitive, all occurrences)
        eval_expr = re.sub(r'\b' + kegg_id + r'\b', str(presence), eval_expr)
    
    # Replace boolean operators with Python equivalents
    eval_expr = eval_expr.replace(' AND ', ' and ')
    eval_expr = eval_expr.replace(' OR ', ' or ')
    eval_expr = eval_expr.replace(' NOT ', ' not ')
    
    # Evaluate the expression
    try:
        result = eval(eval_expr)
        return bool(result)
    except Exception as e:
        print(f"Error evaluating expression '{expression}' -> '{eval_expr}': {e}")
        return False


def evaluate_pathways(
    matrix: pd.DataFrame,
    pathway_defs: pd.DataFrame,
    pathway_genes: dict[str, Set[str]]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate pathway presence and completeness for all genomes.
    
    Args:
        matrix: KEGG gene x genome matrix
        pathway_defs: Pathway definitions with boolean expressions
        pathway_genes: Dictionary mapping pathway to set of gene IDs
    
    Returns:
        Tuple of (presence_df, completeness_df):
          - presence_df: Boolean DataFrame based on marker gene expressions (pathways from pathway_defs)
          - completeness_df: Float DataFrame with fraction of pathway genes present (all pathways)
    """
    genomes = matrix.columns
    
    # Convert matrix index to string for consistent matching
    matrix.index = matrix.index.astype(str)
    
    # ===== PRESENCE CALCULATION =====
    # Only for pathways in pathway_definitions.csv
    pathways = pathway_defs["pathway"].values
    expressions = pathway_defs["expression"].values
    
    print(f"Evaluating presence for {len(pathways)} pathways across {len(genomes)} genomes")
    
    presence_results = []
    
    with tqdm(total=len(pathways) * len(genomes), desc="Evaluating pathway presence") as pbar:
        for pathway, expression in zip(pathways, expressions):
            pathway_presence = {}
            
            for genome in genomes:
                gene_presence = create_gene_presence_dict(matrix, genome)
                pathway_presence[genome] = evaluate_expression(expression, gene_presence)
                pbar.update(1)
            
            presence_results.append({
                "pathway": pathway,
                **pathway_presence
            })
    
    presence_df = pd.DataFrame(presence_results)
    presence_df = presence_df.set_index("pathway")
    
    # ===== COMPLETENESS CALCULATION =====
    # For ALL pathways in kegg_genes.csv
    pathways_list = list(pathway_genes.keys())
    
    print(f"\nCalculating completeness for {len(pathways_list)} pathways across {len(genomes)} genomes")
    
    completeness_results = []
    
    with tqdm(total=len(pathways_list) * len(genomes), desc="Calculating completeness") as pbar:
        for pathway in pathways_list:
            pathway_gene_set = pathway_genes[pathway]
            total_genes = len(pathway_gene_set)
            
            pathway_completeness = {}
            
            for genome in genomes:
                gene_presence = create_gene_presence_dict(matrix, genome)
                
                # Completeness from all pathway genes
                if total_genes > 0:
                    present_count = sum(1 for gene in pathway_gene_set if gene_presence.get(gene, False))
                    pathway_completeness[genome] = present_count / total_genes
                else:
                    pathway_completeness[genome] = 0.0
                
                pbar.update(1)
            
            completeness_results.append({
                "pathway": pathway,
                **pathway_completeness
            })
    
    completeness_df = pd.DataFrame(completeness_results)
    completeness_df = completeness_df.set_index("pathway")
    
    # ===== SUMMARY STATISTICS =====
    genomes_per_pathway = presence_df.sum(axis=1).sort_values(ascending=False)
    
    print(f"\nGenomes per pathway (boolean presence from markers):")
    for pathway, count in genomes_per_pathway.items():
        print(f"  {pathway}: {count}")
    
    avg_completeness = completeness_df.mean(axis=1).sort_values(ascending=False)
    print(f"\nAverage completeness per pathway (from all genes):")
    for pathway, frac in avg_completeness.items():
        print(f"  {pathway}: {frac:.2f}")
    
    return presence_df, completeness_df


def write_output(
    presence_df: pd.DataFrame,
    completeness_df: pd.DataFrame,
    output_path: str,
    output_format: str = "hdf5"
) -> None:
    """Write pathway presence and completeness matrices to file.
    
    Args:
        presence_df: Boolean pathway x genome presence matrix
        completeness_df: Float pathway x genome completeness matrix (0-1)
        output_path: Path to output file
        output_format: Output format ('hdf5' or 'csv')
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True) if os.path.dirname(output_path) else None
    
    if output_format == "hdf5":
        # Convert boolean to int for HDF5 storage
        presence_int = presence_df.astype(int)
        print(f"Writing pathway presence matrix to {output_path}")
        presence_int.to_hdf(output_path, key="pathway_presence", mode="w")
        print(f"Writing pathway completeness matrix to {output_path}")
        completeness_df.to_hdf(output_path, key="pathway_completeness", mode="a")
        
        # Save genome IDs (column names) and pathway names (row names) for coordinate mapping
        with h5py.File(output_path, "a") as f:
            str_dtype = h5py.string_dtype(encoding="utf-8")
            f.create_dataset("genome_ids", data=np.array(presence_df.columns, dtype=str_dtype))
            f.create_dataset("pathway_names", data=np.array(presence_df.index, dtype=str_dtype))
        print(f"Saved genome IDs and pathway names to {output_path}")
    else:
        # For CSV, write separate files
        presence_path = output_path.replace(".csv", "_presence.csv")
        completeness_path = output_path.replace(".csv", "_completeness.csv")
        print(f"Writing pathway presence matrix to {presence_path}")
        presence_df.to_csv(presence_path)
        print(f"Writing pathway completeness matrix to {completeness_path}")
        completeness_df.to_csv(completeness_path)


def build_arg_parser() -> argparse.ArgumentParser:
    """Build argument parser for command-line interface."""
    parser = argparse.ArgumentParser(
        description="Evaluate pathway presence using boolean logic over KEGG genes."
    )
    parser.add_argument(
        "-m", "--matrix",
        required=True,
        help="Input HDF5 matrix file (KEGG genes x genomes)"
    )
    parser.add_argument(
        "-p", "--pathways",
        required=True,
        help="CSV file with pathway definitions (columns: pathway, expression, meaning)"
    )
    parser.add_argument(
        "-g", "--genes",
        required=True,
        help="CSV file with KEGG genes and pathway column for completeness calculation"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file path (HDF5 or CSV based on extension)"
    )
    return parser


def main() -> None:
    """Parse arguments, evaluate pathways, and write output."""
    parser = build_arg_parser()
    args = parser.parse_args()
    
    # Determine output format from extension
    output_format = "hdf5" if args.output.endswith(".h5") else "csv"
    
    matrix = load_matrix(args.matrix)
    pathway_defs = load_pathway_definitions(args.pathways)
    pathway_genes = load_kegg_genes_by_pathway(args.genes)
    
    presence_df, completeness_df = evaluate_pathways(matrix, pathway_defs, pathway_genes)
    
    write_output(presence_df, completeness_df, args.output, output_format)


if __name__ == "__main__":
    main()
