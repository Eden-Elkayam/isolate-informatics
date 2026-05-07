"""
kofamscan_to_pipeline_input.py

Converts KofamScan output TSV files to the .tsv.gz format expected by
c1_informatics (collate_globDB_kegg_annotations.py).

Takes only significant hits (lines marked with * in the first column).

Output: one <genome_id>.tsv.gz per genome, with columns: source, accession

Usage:
    python kofamscan_to_pipeline_input.py \
        --input-dir runs/run_1/results/kofamscan/ \
        --output-dir runs/run_1/results/c1_input/
"""

import argparse
import gzip
import csv
from pathlib import Path


def parse_kofamscan(tsv_path):
    """
    Parse a KofamScan detail-tsv output file.
    Returns list of KO accessions for significant hits only (marked with *).
    """
    rows = []
    with open(tsv_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            significance = parts[0].strip()
            if significance != '*':
                continue
            ko = parts[2].strip()
            if ko.startswith('K'):
                rows.append({
                    "source": "KOfam",
                    "accession": ko
                })
    return rows


def write_tsv_gz(rows, out_path):
    with gzip.open(out_path, 'wt') as f:
        writer = csv.DictWriter(f, fieldnames=["source", "accession"],
                                delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Convert KofamScan TSV output to c1_informatics input format."
    )
    parser.add_argument("--input-dir", required=True,
                        help="Directory containing KofamScan .ko.tsv files")
    parser.add_argument("--output-dir", required=True,
                        help="Directory to write .tsv.gz output files")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob("*.ko.tsv"))
    if not files:
        print(f"No .ko.tsv files found in {input_dir}")
        return

    print(f"Found {len(files)} KofamScan result files\n")

    for tsv_file in files:
        # Derive genome ID by stripping .ko.tsv
        genome_id = tsv_file.name.replace(".ko.tsv", "")
        rows = parse_kofamscan(tsv_file)
        out_path = output_dir / f"{genome_id}.tsv.gz"
        write_tsv_gz(rows, out_path)
        print(f"  {genome_id}: {len(rows)} significant KO hits -> {out_path.name}")

    print(f"\nDone. Output files in: {output_dir}")


if __name__ == "__main__":
    main()
