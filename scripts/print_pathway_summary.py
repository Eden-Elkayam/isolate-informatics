#!/usr/bin/env python3
"""
Print pathway presence matrix: genomes as rows, pathways as columns.
Values are 1 (present) or 0 (absent).

Usage:
  python scripts/print_pathway_summary.py
  python scripts/print_pathway_summary.py --save   # also saves to output/pathway_summary.xlsx
"""

import argparse
import h5py
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--save", action="store_true", help="Save to Excel")
parser.add_argument("--pivot", action="store_true", help="Pivot: pathways as rows, genomes as columns")
parser.add_argument("--output-dir", default="output", help="Output directory (default: output)")
parser.add_argument("--taxonomy", default=None, help="Path to taxonomy_summary.csv; adds a species row")
args = parser.parse_args()

PATHWAY_MATRIX = f"{args.output_dir}/pathway_presence_matrix.h5"

with h5py.File(PATHWAY_MATRIX, "r") as f:
    genomes  = [g.decode() if isinstance(g, bytes) else g for g in f["genome_ids"][:]]
    pathways = [p.decode() if isinstance(p, bytes) else p for p in f["pathway_names"][:]]
    presence = f["pathway_presence/block0_values"][:]

df = pd.DataFrame(presence, index=pathways, columns=genomes).T  # genomes as rows
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

print("\n=== PATHWAY PRESENCE (1=present, 0=absent) ===\n")
print(out_df.to_string())
print()

if args.save:
    out = f"{args.output_dir}/pathway_summary.xlsx"
    out_df.to_excel(out)
    print(f"Saved to {out}")
