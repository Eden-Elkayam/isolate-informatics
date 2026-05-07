#!/usr/bin/env python3
"""
Print KOs found per pathway per genome: genomes as rows, pathways as columns.
Values are comma-separated KO IDs found, or "None" if none found.
Note: shows KOs found regardless of whether the full pathway was detected.

Usage:
  python scripts/print_genes_summary.py
  python scripts/print_genes_summary.py --save
  python scripts/print_genes_summary.py --pivot
  python scripts/print_genes_summary.py --output-dir runs/run_2/output
"""

import re
import argparse
import h5py
import pandas as pd

PATHWAY_DEFS = "data/pathway_definitions.csv"

parser = argparse.ArgumentParser()
parser.add_argument("--save", action="store_true", help="Save to Excel")
parser.add_argument("--pivot", action="store_true", help="Pivot: pathways as rows, genomes as columns")
parser.add_argument("--output-dir", default="output", help="Output directory (default: output)")
parser.add_argument("--taxonomy", default=None, help="Path to taxonomy_summary.csv; adds a species row")
args = parser.parse_args()

PATHWAY_MATRIX = f"{args.output_dir}/pathway_presence_matrix.h5"
GENE_MATRIX    = f"{args.output_dir}/kegg_by_genome_matrix.h5"

with h5py.File(PATHWAY_MATRIX, "r") as f:
    genomes  = [g.decode() if isinstance(g, bytes) else g for g in f["genome_ids"][:]]
    pathways = [p.decode() if isinstance(p, bytes) else p for p in f["pathway_names"][:]]

gene_df = pd.read_hdf(GENE_MATRIX, key="kegg_matrix")

defs = pd.read_csv(PATHWAY_DEFS)
pathway_col = defs.columns[0]
expr_col    = defs.columns[2]

pathway_genes = {}
for _, row in defs.iterrows():
    kos = re.findall(r"K\d{5}", str(row[expr_col]))
    pathway_genes[row[pathway_col]] = kos

records = []
for genome in genomes:
    row = {"genome": genome}
    for pathway in pathways:
        kos = pathway_genes.get(pathway, [])
        found = [ko for ko in kos if ko in gene_df.index and gene_df.loc[ko, genome] > 0]
        row[pathway] = ", ".join(found) if found else "None"
    records.append(row)

df = pd.DataFrame(records).set_index("genome")
out_df = df.T if args.pivot else df

if args.taxonomy:
    tax = pd.read_csv(args.taxonomy).set_index("genome")["species"]
    if args.pivot:
        species_row = pd.DataFrame(
            {g: tax.get(g, "unknown") for g in out_df.columns},
            index=["species"]
        )
        out_df = pd.concat([species_row, out_df])
    else:
        out_df.insert(0, "species", [tax.get(g, "unknown") for g in out_df.index])

print("\n=== KOs FOUND PER PATHWAY ===\n")
print(out_df.to_string())
print()

if args.save:
    out = f"{args.output_dir}/genes_summary.xlsx"
    out_df.to_excel(out)
    print(f"Saved to {out}")
