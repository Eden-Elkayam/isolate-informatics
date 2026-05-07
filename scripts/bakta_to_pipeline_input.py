#!/usr/bin/env python3
"""
bakta_to_pipeline_input.py

Extracts KEGG KO annotations from merged Bakta .gb files and writes
.tsv.gz files in the format expected by c1_informatics.

Does not use BioPython -- uses the same simple line parser as the
rest of the kegg_preprocessing pipeline.

Usage:
    python scripts/bakta_to_pipeline_input.py \
        --input-dir results/merged/ \
        --output-dir results/bakta_input/
"""

import re
import gzip
import argparse
import csv
from pathlib import Path


def extract_kegg_from_gb(gb_path):
    rows = []
    with open(gb_path) as f:
        for line in f:
            line_stripped = line.strip()
            if line_stripped.startswith('/db_xref="KEGG:'):
                ko = re.search(r'KEGG:(K\d+)', line_stripped)
                if ko:
                    rows.append({
                        "source": "KOfam",
                        "accession": ko.group(1)
                    })
    return rows


def write_tsv_gz(rows, out_path):
    with gzip.open(out_path, 'wt') as f:
        writer = csv.DictWriter(f, fieldnames=["source", "accession"], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Extract Bakta KEGG annotations from merged .gb files."
    )
    parser.add_argument("--input-dir", required=True,
                        help="Directory containing merged .gb files")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for .tsv.gz files")
    parser.add_argument("--ext", default="gb",
                        help="File extension (default: gb)")
    args = parser.parse_args()

    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob(f"*.{args.ext}"))
    if not files:
        print(f"No .{args.ext} files found in {input_dir}")
        return

    print(f"Found {len(files)} genome(s)\n")

    for gb_file in files:
        genome_id = gb_file.stem
        rows = extract_kegg_from_gb(gb_file)
        out_path = output_dir / f"{genome_id}.tsv.gz"
        write_tsv_gz(rows, out_path)
        print(f"  {genome_id}: {len(rows)} KEGG annotations -> {out_path.name}")

    print(f"\nDone. Files written to {output_dir}")


if __name__ == "__main__":
    main()
