#!/usr/bin/env python3
"""
export_excel.py

Exports a 3-sheet Excel summary for a completed pipeline run:
  1. Phylogeny  — full GTDB lineage per genome
  2. Pathways   — pathway presence (0/1) per genome, with phylogeny columns
  3. Genes      — KOs found per pathway per genome, with phylogeny columns

Usage:
    python scripts/export_excel.py \
        --output-dir runs/<run>/pathways \
        --taxonomy   runs/<run>/phylogeny/taxonomy/taxonomy_summary.csv \
        --output     runs/<run>/summary.xlsx
"""

import re
import argparse
import h5py
import pandas as pd
from pathlib import Path

PATHWAY_DEFS = "data/pathway_definitions.csv"
PHYLO_COLS   = ["species", "genus", "family", "full_lineage"]


def load_taxonomy(taxonomy_path):
    """Load taxonomy CSV; return DataFrame indexed by genome."""
    tax = pd.read_csv(taxonomy_path)
    # Ensure all expected columns are present
    for col in PHYLO_COLS:
        if col not in tax.columns:
            tax[col] = "unknown"
    return tax.set_index("genome")


def load_pathway_matrix(output_dir):
    """Return (genomes list, pathways list, presence ndarray)."""
    path = Path(output_dir) / "pathway_presence_matrix.h5"
    with h5py.File(path, "r") as f:
        genomes  = [g.decode() if isinstance(g, bytes) else g for g in f["genome_ids"][:]]
        pathways = [p.decode() if isinstance(p, bytes) else p for p in f["pathway_names"][:]]
        presence = f["pathway_presence/block0_values"][:]
    return genomes, pathways, presence


def load_gene_matrix(output_dir):
    """Return KEGG x genome DataFrame."""
    path = Path(output_dir) / "kegg_by_genome_matrix.h5"
    return pd.read_hdf(path, key="kegg_matrix")


def build_phylo_prefix(genomes, tax_df):
    """Build a DataFrame with phylogeny columns for the given genome list."""
    rows = []
    for g in genomes:
        if g in tax_df.index:
            row = tax_df.loc[g, PHYLO_COLS].to_dict()
        else:
            row = {col: "unknown" for col in PHYLO_COLS}
        rows.append(row)
    return pd.DataFrame(rows, index=genomes)


def build_pathways_df(genomes, pathways, presence, tax_df):
    pathway_df = pd.DataFrame(presence, index=pathways, columns=genomes).T
    phylo = build_phylo_prefix(genomes, tax_df)
    return pd.concat([phylo, pathway_df], axis=1)


def write_genes_sheet(writer, genomes, pathways, gene_df, tax_df):
    """
    Write the Genes sheet directly via openpyxl with explicit two-row headers:
      Row 1 : "genome", phylo col names, then pathway name repeated once per KO
      Row 2 : empty for genome/phylo cols, KO ID for each data column
      Row 3+: genome_id, phylo values, True/False per KO
    """
    defs = pd.read_csv(PATHWAY_DEFS)
    pathway_col = defs.columns[0]
    expr_col    = defs.columns[2]

    pathway_genes = {}
    for _, row in defs.iterrows():
        kos = re.findall(r"K\d{5}", str(row[expr_col]))
        pathway_genes[row[pathway_col]] = kos

    ko_columns = []  # [(pathway, ko), ...]
    for pathway in pathways:
        for ko in pathway_genes.get(pathway, []):
            ko_columns.append((pathway, ko))

    META = 1 + len(PHYLO_COLS)  # genome col + phylo cols

    ws = writer.book.create_sheet("Genes")

    # Row 1: genome, phylo names, pathway name repeated per KO column
    ws.cell(row=1, column=1, value="genome")
    for i, col in enumerate(PHYLO_COLS):
        ws.cell(row=1, column=2 + i, value=col)
    for i, (pathway, _) in enumerate(ko_columns):
        ws.cell(row=1, column=META + 1 + i, value=pathway)

    # Row 2: KO IDs (genome/phylo cols left blank)
    for i, (_, ko) in enumerate(ko_columns):
        ws.cell(row=2, column=META + 1 + i, value=ko)

    # Data rows (row 3+)
    for r, genome in enumerate(genomes):
        if genome in tax_df.index:
            phylo_vals = [tax_df.loc[genome, col] for col in PHYLO_COLS]
        else:
            phylo_vals = ["unknown"] * len(PHYLO_COLS)

        ws.cell(row=3 + r, column=1, value=genome)
        for i, val in enumerate(phylo_vals):
            ws.cell(row=3 + r, column=2 + i, value=val)
        for i, (pathway, ko) in enumerate(ko_columns):
            present = (
                ko in gene_df.index
                and genome in gene_df.columns
                and gene_df.loc[ko, genome] > 0
            )
            ws.cell(row=3 + r, column=META + 1 + i, value=int(present))

    return ko_columns


def main():
    parser = argparse.ArgumentParser(description="Export 3-sheet Excel summary.")
    parser.add_argument("--output-dir", required=True,
                        help="Pathways output directory (contains pathway_presence_matrix.h5)")
    parser.add_argument("--taxonomy", default=None,
                        help="Path to taxonomy_summary.csv (optional)")
    parser.add_argument("--output", required=True,
                        help="Output .xlsx file path")
    args = parser.parse_args()

    genomes, pathways, presence = load_pathway_matrix(args.output_dir)
    gene_df = load_gene_matrix(args.output_dir)

    # Load taxonomy or fall back to unknowns
    if args.taxonomy and Path(args.taxonomy).exists():
        tax_df = load_taxonomy(args.taxonomy)
    else:
        tax_df = pd.DataFrame(
            {col: ["unknown"] * len(genomes) for col in PHYLO_COLS},
            index=genomes
        )

    # Sheet 1: Phylogeny
    phylo_sheet = tax_df.loc[
        [g for g in genomes if g in tax_df.index]
    ][PHYLO_COLS].copy()
    phylo_sheet.index.name = "genome"
    phylo_sheet = phylo_sheet.reset_index()

    # Sheet 2: Pathways (genomes as rows)
    pathways_sheet = build_pathways_df(genomes, pathways, presence, tax_df)
    pathways_sheet.index.name = "genome"
    pathways_sheet = pathways_sheet.reset_index()

    # Write Excel — Phylogeny and Pathways via pandas; Genes via openpyxl directly
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        phylo_sheet.to_excel(writer,    sheet_name="Phylogeny", index=False)
        pathways_sheet.to_excel(writer, sheet_name="Pathways",  index=False)
        ko_columns = write_genes_sheet(writer, genomes, pathways, gene_df, tax_df)

    print(f"Saved: {args.output}")
    print(f"  Sheet 1 — Phylogeny : {len(phylo_sheet)} genome(s)")
    print(f"  Sheet 2 — Pathways  : {len(pathways_sheet)} genome(s) x {len(pathways)} pathways")
    print(f"  Sheet 3 — Genes     : {len(genomes)} genome(s) x {len(ko_columns)} KO columns ({len(pathways)} pathways)")


if __name__ == "__main__":
    main()
