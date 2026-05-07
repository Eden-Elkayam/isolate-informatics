#!/usr/bin/env python3
"""
sourmash_taxonomy.py

Runs sourmash sketch + gather + tax genome on each .fna genome file to
identify species using the GTDB database with proper GTDB taxonomy lineages.

Outputs:
  - results/taxonomy/gtdb_taxonomy.lineage.csv  — full GTDB lineage per genome
  - results/taxonomy/taxonomy_summary.csv        — simplified species name per genome
  - results/taxonomy/<genome>.gather.csv         — full sourmash gather output

Usage:
    python scripts/sourmash_taxonomy.py \
        --fna-dir results/fna/ \
        --db ~/db/sourmash/gtdb-rs214-reps.k31.zip \
        --lineages ~/db/sourmash/taxonomy/gtdb-rs214.lineages.csv \
        --output-dir results/taxonomy/
"""

import argparse
import subprocess
import csv
from pathlib import Path


def run_cmd(cmd, description=""):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {description}")
        print(f"  {result.stderr.strip()}")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Run sourmash taxonomy on genome FNA files using GTDB."
    )
    parser.add_argument("--fna-dir", required=True,
                        help="Directory containing .fna genome files")
    parser.add_argument("--db", required=True,
                        help="Sourmash GTDB database zip file")
    parser.add_argument("--lineages", required=True,
                        help="GTDB lineages CSV for sourmash tax genome "
                             "(e.g. gtdb-rs214.lineages.csv)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for taxonomy results")
    args = parser.parse_args()

    fna_dir    = Path(args.fna_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fna_files = sorted(fna_dir.glob("*.fna"))
    if not fna_files:
        print(f"No .fna files found in {fna_dir}")
        return

    print(f"Found {len(fna_files)} genome(s)\n")

    gather_csvs = []

    for fna_file in fna_files:
        genome_id  = fna_file.stem
        sig_file   = output_dir / f"{genome_id}.sig"
        gather_csv = output_dir / f"{genome_id}.gather.csv"

        print(f"  Processing {genome_id} ...")

        # Step 1: sketch
        ok = run_cmd(
            f"sourmash sketch dna -p k=31,scaled=1000 {fna_file} "
            f"-o {sig_file} --name {genome_id}",
            "sourmash sketch"
        )
        if not ok:
            print(f"  Skipping {genome_id} — sketch failed")
            continue

        # Step 2: gather
        ok = run_cmd(
            f"sourmash gather {sig_file} {args.db} "
            f"-o {gather_csv} --threshold-bp 50000 -q",
            "sourmash gather"
        )
        if not ok or not gather_csv.exists():
            print(f"  No match found for {genome_id}")
            continue

        gather_csvs.append(str(gather_csv))
        print(f"  Gather complete: {gather_csv.name}")

    if not gather_csvs:
        print("No gather results to classify.")
        return

    # Step 3: sourmash tax genome — assign GTDB lineages
    print(f"\nRunning sourmash tax genome on {len(gather_csvs)} genome(s) ...")
    gather_args = " \\\n    -g ".join(gather_csvs)
    lineage_out = output_dir / "gtdb_taxonomy"

    ok = run_cmd(
        f"sourmash tax genome \\\n    -g {gather_args} \\\n"
        f"    -t {args.lineages} \\\n"
        f"    -o {lineage_out} \\\n"
        f"    -F lineage_csv",
        "sourmash tax genome"
    )
    if not ok:
        print("Tax genome step failed.")
        return

    lineage_csv = Path(str(lineage_out) + ".lineage.csv")

    # Step 4: parse lineage CSV into simple summary
    summary_path = output_dir / "taxonomy_summary.csv"
    results = []
    if lineage_csv.exists():
        with open(lineage_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                results.append({
                    "genome": row["ident"],
                    "species": row.get("species", "unknown"),
                    "genus": row.get("genus", "unknown"),
                    "family": row.get("family", "unknown"),
                    "full_lineage": ";".join([
                        row.get("superkingdom", ""),
                        row.get("phylum", ""),
                        row.get("class", ""),
                        row.get("order", ""),
                        row.get("family", ""),
                        row.get("genus", ""),
                        row.get("species", "")
                    ])
                })

    with open(summary_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["genome", "species", "genus",
                                                "family", "full_lineage"])
        writer.writeheader()
        writer.writerows(results)

    print(f"\nDone. Results saved to: {output_dir}")
    print(f"\n=== GTDB TAXONOMY SUMMARY ===\n")
    for r in results:
        print(f"  {r['genome']:<25} {r['species']}")
    print()


if __name__ == "__main__":
    main()
