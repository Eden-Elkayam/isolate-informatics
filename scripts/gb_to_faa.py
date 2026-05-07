"""
gb_to_faa.py

Converts merged .gb annotation files to .faa protein FASTA files
ready for input to KofamScan.

Usage:
    python gb_to_faa.py <input_dir> [--out <output_dir>] [--ext gb]
"""

import re
import argparse
from pathlib import Path


def extract_proteins(gb_path):
    """
    Parses a GenBank file and extracts locus_tag + translation for every CDS.
    Prefixes duplicate locus tags with a counter to ensure unique names.
    Returns a list of (header, sequence) tuples.
    """
    proteins = []
    seen = {}
    current = None
    in_translation = False
    seq_buffer = []

    def unique_tag(tag):
        if tag not in seen:
            seen[tag] = 0
            return tag
        seen[tag] += 1
        return f"{tag}_{seen[tag]}"

    with open(gb_path) as f:
        for line in f:
            line_stripped = line.strip()

            # Detect CDS feature
            if re.match(r'^\s{5}CDS\s+', line):
                # Save previous if complete
                if current and seq_buffer:
                    proteins.append((unique_tag(current), ''.join(seq_buffer).rstrip('*')))
                current = None
                in_translation = False
                seq_buffer = []
                continue

            # Pick up locus_tag
            if line_stripped.startswith('/locus_tag='):
                tag = re.sub(r'^/locus_tag="?|"?$', '', line_stripped)
                current = tag
                continue

            # Pick up label as fallback if no locus_tag
            if line_stripped.startswith('/label=') and not current:
                label = re.sub(r'^/label="?|"?$', '', line_stripped)
                current = label
                continue

            # Start of translation
            if line_stripped.startswith('/translation='):
                in_translation = True
                seq = re.sub(r'^/translation="?', '', line_stripped).rstrip('"')
                seq_buffer.append(seq)
                if line_stripped.endswith('"'):
                    in_translation = False
                continue

            # Continuation lines of translation
            if in_translation:
                seq_buffer.append(line_stripped.rstrip('"'))
                if line_stripped.endswith('"'):
                    in_translation = False
                continue

        # Don't forget the last CDS in the file
        if current and seq_buffer:
            proteins.append((unique_tag(current), ''.join(seq_buffer).rstrip('*')))

    return proteins


def write_faa(proteins, out_path):
    with open(out_path, 'w') as f:
        for header, seq in proteins:
            f.write(f">{header}\n")
            # Wrap at 60 chars
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def main():
    parser = argparse.ArgumentParser(description="Convert merged GB files to FAA for KofamScan.")
    parser.add_argument("input_dir", help="Directory containing merged .gb files")
    parser.add_argument("--out", help="Output directory for .faa files (default: same as input)", default=None)
    parser.add_argument("--ext", default="gb", help="File extension to look for (default: gb)")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    out_dir = Path(args.out) if args.out else input_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob(f"*.{args.ext}"))
    if not files:
        print(f"No .{args.ext} files found in {input_dir}")
        return

    print(f"Found {len(files)} .{args.ext} files\n")

    for gb_file in files:
        proteins = extract_proteins(gb_file)
        out_path = out_dir / (gb_file.stem + ".faa")
        write_faa(proteins, out_path)
        print(f"  {gb_file.name} -> {out_path.name}  ({len(proteins)} proteins)")

    print(f"\nDone. FAA files written to {out_dir}")
    print("\nNext step — run KofamScan on each .faa file:")
    print(f"  exec_annotation -p profiles/ -k ko_list --cpu 8 -f detail-tsv -o <genome>.ko.tsv <genome>.faa")


if __name__ == "__main__":
    main()
