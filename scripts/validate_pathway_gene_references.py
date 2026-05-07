#!/usr/bin/env python3
"""Validate KEGG IDs referenced in pathway definitions.

Checks that all KEGG IDs referenced in pathway boolean expressions are present
in the KEGG gene reference table.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


KEGG_PATTERN = re.compile(r"K\d{5}")


def _resolve_default_pathways_file(path: Path) -> Path:
    """Support either pathway_definitions.csv or pathways_definitions.csv."""
    if path.exists():
        return path

    alt = path.parent / "pathways_definitions.csv"
    if alt.exists():
        return alt

    return path


def load_genes(genes_csv: Path) -> set[str]:
    """Load KEGG IDs from kegg_genes.csv."""
    with genes_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {genes_csv}")

        preferred_columns = ["KO", "ko", "kegg_id", "gene"]
        gene_column = next((col for col in preferred_columns if col in reader.fieldnames), None)
        if gene_column is None:
            gene_column = reader.fieldnames[0]

        genes = {
            row[gene_column].strip()
            for row in reader
            if row.get(gene_column) and row[gene_column].strip()
        }

    return genes


def parse_pathway_rows(pathways_csv: Path) -> list[tuple[str, str]]:
    """Return (pathway_name, expression) rows from pathway definitions."""
    with pathways_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {pathways_csv}")

        required = {"pathway", "expression"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Missing required columns in {pathways_csv}: {', '.join(sorted(missing))}"
            )

        rows = []
        for row in reader:
            rows.append((row["pathway"].strip(), (row["expression"] or "").strip()))

    return rows


def extract_kegg_ids(expression: str) -> set[str]:
    """Extract all K##### identifiers from an expression string."""
    return set(KEGG_PATTERN.findall(expression))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate that pathway expression KEGG IDs exist in kegg_genes.csv"
    )
    parser.add_argument(
        "--genes",
        default="data/kegg_genes.csv",
        help="Path to kegg_genes.csv (default: data/kegg_genes.csv)",
    )
    parser.add_argument(
        "--pathways",
        default="data/pathway_definitions.csv",
        help=(
            "Path to pathway definitions CSV (default: data/pathway_definitions.csv; "
            "also auto-checks data/pathways_definitions.csv if default is missing)"
        ),
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    genes_path = Path(args.genes)
    pathways_path = _resolve_default_pathways_file(Path(args.pathways))

    if not genes_path.exists():
        raise FileNotFoundError(f"Genes file not found: {genes_path}")
    if not pathways_path.exists():
        raise FileNotFoundError(f"Pathways file not found: {pathways_path}")

    valid_genes = load_genes(genes_path)
    pathway_rows = parse_pathway_rows(pathways_path)

    missing_by_pathway: dict[str, set[str]] = {}
    all_referenced: set[str] = set()

    for pathway_name, expression in pathway_rows:
        referenced = extract_kegg_ids(expression)
        all_referenced.update(referenced)

        missing = referenced - valid_genes
        if missing:
            missing_by_pathway[pathway_name] = missing

    print(f"Genes in reference file: {len(valid_genes)}")
    print(f"Unique KEGG IDs referenced in expressions: {len(all_referenced)}")

    if not missing_by_pathway:
        print("OK: all KEGG IDs referenced in pathway expressions are present in genes file.")
        return

    print("\nMissing KEGG IDs detected:")
    for pathway_name in sorted(missing_by_pathway):
        missing = ", ".join(sorted(missing_by_pathway[pathway_name]))
        print(f"  - {pathway_name}: {missing}")

    raise SystemExit(1)


if __name__ == "__main__":
    main()
