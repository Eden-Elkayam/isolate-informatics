"""
find_missing_kegg.py

Scans a GenBank-format annotation file for genes that match
formate dehydrogenase (or any target pathway) by product name
but lack a KEGG KO assignment.

Usage:
    python find_missing_kegg.py <annotation.gbk> [--terms terms.txt]
"""

import re
import sys
import argparse
from pathlib import Path

# ── Target terms ─────────────────────────────────────────────────────────────
# These mirror the gene names used in c1_informatics KEGG expressions for FDH.
# Extend this list to cover other pathways of interest.

FDH_TERMS = [
    r"formate dehydrogenase",
    r"fdnG",
    r"fdnH",
    r"fdnI",
    r"fdhF",
    r"fdoG",
    r"fdoH",
    r"fdoI",
    r"formate:acceptor oxidoreductase",
    r"formate oxidoreductase",
]

# ── Parser ────────────────────────────────────────────────────────────────────

def parse_features(gbk_path):
    """
    Minimal GenBank CDS parser.
    Returns a list of dicts, one per CDS, with keys:
        locus_tag, product, label, kegg_ids, has_kegg, coords
    """
    features = []
    current = None

    with open(gbk_path) as f:
        for line in f:
            line_stripped = line.strip()

            # Detect start of a CDS feature
            cds_match = re.match(r'^\s{5}CDS\s+(complement\()?(\d+)\.\.(\d+)', line)
            if cds_match:
                current = {
                    "coords": line_stripped,
                    "locus_tag": "",
                    "product": "",
                    "label": "",
                    "kegg_ids": [],
                    "has_kegg": False,
                }
                features.append(current)
                continue

            if current is None:
                continue

            # Stop collecting if we hit a new feature key
            if re.match(r'^\s{5}\w', line) and not line_stripped.startswith('/'):
                current = None
                continue

            # Parse qualifiers
            if line_stripped.startswith('/product='):
                current["product"] = re.sub(r'^/product="?|"?$', '', line_stripped)

            elif line_stripped.startswith('/locus_tag='):
                current["locus_tag"] = re.sub(r'^/locus_tag="?|"?$', '', line_stripped)

            elif line_stripped.startswith('/label='):
                current["label"] = re.sub(r'^/label="?|"?$', '', line_stripped)

            elif line_stripped.startswith('/db_xref="KEGG:'):
                ko = re.search(r'KEGG:(K\d+)', line_stripped)
                if ko:
                    current["kegg_ids"].append(ko.group(1))
                    current["has_kegg"] = True

    return features


def find_unrepresented(features, terms):
    """
    Returns CDS entries whose product/label matches any term
    but that lack a KEGG KO assignment.
    """
    pattern = re.compile('|'.join(terms), re.IGNORECASE)
    missing = []
    matched = []

    for feat in features:
        search_text = f"{feat['product']} {feat['label']}"
        if pattern.search(search_text):
            matched.append(feat)
            if not feat["has_kegg"]:
                missing.append(feat)

    return matched, missing


def report(matched, missing, terms):
    has_kegg = [f for f in matched if f["has_kegg"]]

    print(f"\n{'='*60}")
    print(f"  Target terms: {', '.join(terms[:3])}{'...' if len(terms) > 3 else ''}")
    print(f"{'='*60}")
    print(f"  Total CDS matching target terms : {len(matched)}")
    print(f"  Of those, WITH KEGG KO          : {len(has_kegg)}")
    print(f"  Of those, MISSING KEGG KO       : {len(missing)}")
    print(f"{'='*60}\n")

    if has_kegg:
        print("── Genes WITH KEGG KO ──────────────────────────────────\n")
        for feat in has_kegg:
            tag = feat["locus_tag"] or feat["label"] or "unknown"
            print(f"  Locus tag : {tag}")
            print(f"  Product   : {feat['product']}")
            print(f"  KEGG KO   : {', '.join(feat['kegg_ids'])}")
            print(f"  Location  : {feat['coords']}")
            print()

    if missing:
        print("── Genes MISSING KEGG KO ───────────────────────────────\n")
        for feat in missing:
            tag = feat["locus_tag"] or feat["label"] or "unknown"
            print(f"  Locus tag : {tag}")
            print(f"  Product   : {feat['product']}")
            print(f"  Location  : {feat['coords']}")
            print()
        print(f"Summary: {len(missing)}/{len(matched)} matched genes have no KEGG KO.")
        print("Consider running KofamScan on these proteins to assign KOs upstream.\n")
    else:
        print("No gaps found — all matched genes have KEGG KO assignments.\n")


# ── Main ──────────────────────────────────────────────────────────────────────

def process_one(gb_path, terms, out_dir=None):
    print(f"\nParsing {gb_path.name} ...")
    features = parse_features(gb_path)
    print(f"Found {len(features)} CDS features total.")
    matched, missing = find_unrepresented(features, terms)

    # Print to stdout
    print(f"\n  Genome: {gb_path.stem}")
    report(matched, missing, terms)

    # Optionally write per-genome report file
    if out_dir:
        out_path = out_dir / f"{gb_path.stem}.fdh_report.txt"
        with open(out_path, 'w') as f:
            has_kegg = [feat for feat in matched if feat["has_kegg"]]
            f.write(f"Genome: {gb_path.stem}\n")
            f.write(f"Total FDH-matched CDS : {len(matched)}\n")
            f.write(f"With KEGG KO          : {len(has_kegg)}\n")
            f.write(f"Missing KEGG KO       : {len(missing)}\n\n")

            if has_kegg:
                f.write("── WITH KEGG KO ────────────────────────────\n\n")
                for feat in has_kegg:
                    tag = feat["locus_tag"] or feat["label"] or "unknown"
                    f.write(f"  {tag} | {feat['product']} | {', '.join(feat['kegg_ids'])}\n")
                f.write("\n")

            if missing:
                f.write("── MISSING KEGG KO ─────────────────────────\n\n")
                for feat in missing:
                    tag = feat["locus_tag"] or feat["label"] or "unknown"
                    f.write(f"  {tag} | {feat['product']}\n")

        print(f"  Report saved: {out_path.name}")

    return len(matched), len(missing)


def main():
    parser = argparse.ArgumentParser(description="Find pathway genes missing KEGG KO in Bakta GB output.")
    parser.add_argument("input", help="Path to a single .gb file OR a directory of merged .gb files")
    parser.add_argument("--terms", help="Optional: path to a text file with one regex term per line")
    parser.add_argument("--ext", default="gb", help="File extension to glob when input is a directory (default: gb)")
    parser.add_argument("--out", help="Optional: directory to write per-genome report files")
    args = parser.parse_args()

    terms = FDH_TERMS
    if args.terms:
        with open(args.terms) as f:
            terms = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(terms)} custom terms from {args.terms}")

    out_dir = Path(args.out) if args.out else None
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

    input_path = Path(args.input)

    # Single file mode
    if input_path.is_file():
        process_one(input_path, terms, out_dir)

    # Batch mode — directory of merged files
    elif input_path.is_dir():
        files = sorted(input_path.glob(f"*.{args.ext}"))
        if not files:
            print(f"No .{args.ext} files found in {input_path}")
            return

        print(f"Found {len(files)} .{args.ext} files in {input_path}")
        summary = []

        for gb_file in files:
            n_matched, n_missing = process_one(gb_file, terms, out_dir)
            summary.append((gb_file.stem, n_matched, n_missing))

        # Print summary table
        print(f"\n{'='*60}")
        print(f"  BATCH SUMMARY")
        print(f"{'='*60}")
        print(f"  {'Genome':<30} {'Matched':>8} {'Missing':>8}")
        print(f"  {'-'*46}")
        for name, matched, missing in summary:
            print(f"  {name:<30} {matched:>8} {missing:>8}")
        print(f"{'='*60}\n")

    else:
        print(f"Input not found: {input_path}")


if __name__ == "__main__":
    main()
