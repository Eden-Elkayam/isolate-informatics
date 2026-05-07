#!/usr/bin/env python3
"""
gb_to_fna.py

Extracts nucleotide sequences from merged Bakta .gb files into .fna FASTA
format for use with Sourmash taxonomy classification.

Each contig becomes one FASTA record. The sequence ID is the genome name
plus a contig counter.

Usage:
    python scripts/gb_to_fna.py <input_dir> --out <output_dir> [--ext gb]
"""

import re
import argparse
from pathlib import Path


def extract_sequences(gb_path):
    """
    Parse a merged GB file and extract nucleotide sequences.
    Returns list of (header, sequence) tuples.
    """
    sequences = []
    in_origin = False
    current_seq = []
    contig_count = 0
    genome_name = gb_path.stem

    with open(gb_path) as f:
        for line in f:
            # New record starts with LOCUS
            if line.startswith('LOCUS'):
                if current_seq:
                    seq = ''.join(current_seq)
                    sequences.append((f"{genome_name}_contig{contig_count}", seq))
                    current_seq = []
                contig_count += 1
                in_origin = False

            elif line.startswith('ORIGIN'):
                in_origin = True

            elif line.startswith('//'):
                if current_seq:
                    seq = ''.join(current_seq)
                    sequences.append((f"{genome_name}_contig{contig_count}", seq))
                    current_seq = []
                in_origin = False

            elif in_origin:
                # Strip line numbers and spaces, keep only sequence
                seq_part = re.sub(r'[\d\s]', '', line.strip())
                current_seq.append(seq_part.upper())

    return sequences


def write_fna(sequences, out_path):
    with open(out_path, 'w') as f:
        for header, seq in sequences:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract nucleotide FASTA from merged GB files."
    )
    parser.add_argument("input_dir", help="Directory containing merged .gb files")
    parser.add_argument("--out", required=True, help="Output directory for .fna files")
    parser.add_argument("--ext", default="gb", help="File extension (default: gb)")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob(f"*.{args.ext}"))
    if not files:
        print(f"No .{args.ext} files found in {input_dir}")
        return

    print(f"Found {len(files)} .{args.ext} files\n")

    for gb_file in files:
        sequences = extract_sequences(gb_file)
        out_path = output_dir / (gb_file.stem + ".fna")
        write_fna(sequences, out_path)
        total_bp = sum(len(s) for _, s in sequences)
        print(f"  {gb_file.name} -> {out_path.name}  "
              f"({len(sequences)} contigs, {total_bp:,} bp)")

    print(f"\nDone. FNA files written to {output_dir}")


if __name__ == "__main__":
    main()
