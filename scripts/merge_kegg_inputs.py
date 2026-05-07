#!/usr/bin/env python3
"""
merge_kegg_inputs.py

Merges Bakta and KofamScan KEGG annotations per genome into a single
deduplicated .tsv.gz file for input to c1_informatics.

For each genome, takes the union of KO accessions from both sources.

Usage:
    python scripts/merge_kegg_inputs.py \
        --bakta-dir results/bakta_input/ \
        --kofamscan-dir results/kofamscan_input/ \
        --output-dir results/c1_input/
"""

import argparse
import gzip
import pandas as pd
from pathlib import Path


def load_tsv_gz(path):
    with gzip.open(path, 'rt') as f:
        return pd.read_csv(f, sep='\t')


def write_tsv_gz(df, path):
    with gzip.open(path, 'wt') as f:
        df.to_csv(f, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Merge Bakta and KofamScan KEGG inputs for c1_informatics."
    )
    parser.add_argument("--bakta-dir", required=True,
                        help="Directory containing Bakta .tsv.gz files")
    parser.add_argument("--kofamscan-dir", required=True,
                        help="Directory containing KofamScan .tsv.gz files")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for merged .tsv.gz files")
    args = parser.parse_args()

    bakta_dir     = Path(args.bakta_dir)
    kofamscan_dir = Path(args.kofamscan_dir)
    output_dir    = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Collect all genome IDs from both directories
    bakta_files     = {f.stem.replace('.tsv', ''): f for f in bakta_dir.glob("*.tsv.gz")}
    kofamscan_files = {f.stem.replace('.tsv', ''): f for f in kofamscan_dir.glob("*.tsv.gz")}

    all_genomes = sorted(set(bakta_files) | set(kofamscan_files))

    if not all_genomes:
        print("No .tsv.gz files found in either directory.")
        return

    print(f"Found {len(all_genomes)} genome(s)\n")

    for genome in all_genomes:
        frames = []

        if genome in bakta_files:
            df = load_tsv_gz(bakta_files[genome])
            frames.append(df)
            bakta_count = len(df)
        else:
            bakta_count = 0
            print(f"  WARNING: {genome} missing from Bakta input")

        if genome in kofamscan_files:
            df = load_tsv_gz(kofamscan_files[genome])
            frames.append(df)
            kofamscan_count = len(df)
        else:
            kofamscan_count = 0
            print(f"  WARNING: {genome} missing from KofamScan input")

        merged = pd.concat(frames, ignore_index=True).drop_duplicates()
        out_path = output_dir / f"{genome}.tsv.gz"
        write_tsv_gz(merged, out_path)

        print(f"  {genome}: {bakta_count} Bakta + {kofamscan_count} KofamScan "
              f"-> {len(merged)} unique KOs")

    print(f"\nDone. Merged files in: {output_dir}")


if __name__ == "__main__":
    main()
