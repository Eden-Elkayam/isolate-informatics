#!/usr/bin/env python3
"""
Convert Bakta-annotated GenBank files to the .tsv.gz format expected by
collate_globDB_kegg_annotations.py.

For each .gbk file, extracts KEGG IDs from db_xref fields and writes a
gzipped TSV with columns: source, accession

Output file per genome: <input_dir>/<genome_id>.tsv.gz

Usage:
  python bakta_gbk_to_pipeline_input.py \
    --input-dir /path/to/gbk_files \
    --output-dir /path/to/output

  # For a single file:
  python bakta_gbk_to_pipeline_input.py \
    --input-dir /path/to/gbk_files \
    --output-dir /path/to/output \
    --genome-id my_genome_name
"""

import argparse
import gzip
import os
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def extract_kegg_from_gbk(gbk_path: str) -> list[dict]:
    """Extract KEGG annotations from a Bakta GenBank file.
    
    Returns list of dicts with keys: source, accession
    """
    rows = []
    for record in SeqIO.parse(gbk_path, "genbank"):
        for feat in record.features:
            if feat.type != "CDS":
                continue
            db_xrefs = feat.qualifiers.get("db_xref", [])
            for xref in db_xrefs:
                if xref.startswith("KEGG:"):
                    kegg_id = xref.replace("KEGG:", "").strip()
                    rows.append({
                        "source": "KOfam",
                        "accession": kegg_id,
                    })
    return rows


def genome_id_from_path(gbk_path: str) -> str:
    """Derive genome ID from filename by stripping extension."""
    name = Path(gbk_path).stem
    # Strip common suffixes
    for suffix in ["_annotations", "_bakta", "_genomic"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def convert_gbk_to_tsv_gz(gbk_path: str, output_dir: str, genome_id: str = None) -> str:
    """Convert a single .gbk file to .tsv.gz and return output path."""
    if genome_id is None:
        genome_id = genome_id_from_path(gbk_path)

    rows = extract_kegg_from_gbk(gbk_path)

    if not rows:
        print(f"  WARNING: no KEGG annotations found in {gbk_path}")

    df = pd.DataFrame(rows, columns=["source", "accession"])

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{genome_id}.tsv.gz")

    with gzip.open(out_path, "wt") as f:
        df.to_csv(f, sep="\t", index=False)

    print(f"  {genome_id}: {len(rows)} KEGG annotations → {out_path}")
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Convert Bakta .gbk files to pipeline .tsv.gz input")
    parser.add_argument("--input-dir", required=True, help="Directory containing genome subfolders")
    parser.add_argument("--output-dir", required=True, help="Output directory for .tsv.gz files")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    # Group .gb/.gbk files by their parent folder (= one genome per folder)
    from collections import defaultdict
    genome_files = defaultdict(list)
    for ext in ["*.gb", "*.gbk"]:
        for f in input_dir.glob(f"**/{ext}"):
            # Use the immediate parent folder name as genome ID
            # but only if it's a subfolder of input_dir, not input_dir itself
            if f.parent == input_dir:
                genome_id = f.stem  # fallback: use filename
            else:
                genome_id = f.parent.name
            genome_files[genome_id].append(f)

    if not genome_files:
        print(f"No .gbk or .gb files found in {input_dir}")
        return

    print(f"Found {len(genome_files)} genome(s)")

    for genome_id, files in sorted(genome_files.items()):
        print(f"\nProcessing genome: {genome_id} ({len(files)} contig file(s))")
        all_rows = []
        for f in sorted(files):
            rows = extract_kegg_from_gbk(str(f))
            all_rows.extend(rows)
            print(f"  {f.name}: {len(rows)} KEGG annotations")

        df = pd.DataFrame(all_rows, columns=["source", "accession"]) if all_rows else \
             pd.DataFrame(columns=["source", "accession"])

        out_path = os.path.join(args.output_dir, f"{genome_id}.tsv.gz")
        with gzip.open(out_path, "wt") as fh:
            df.to_csv(fh, sep="\t", index=False)

        print(f"  → {genome_id}: {len(all_rows)} total KEGG annotations → {out_path}")

    print(f"\nDone. Files written to {args.output_dir}")
    print("Next step: set input_dir in config.yaml to this output directory and run the pipeline.")


if __name__ == "__main__":
    main()