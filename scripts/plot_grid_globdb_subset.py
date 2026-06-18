"""
plot_grid_globdb_subset.py

Plot the KO grid for the first run (runs/06032026_5PM/summary.xlsx),
keeping only genomes whose reference was found in GlobDB
(as recorded in runs/globdb_run/representitives.csv).

Usage
-----
python scripts/plot_grid_globdb_subset.py \
    --representatives runs/globdb_run/representitives.csv \
    --summary         runs/06032026_5PM/summary.xlsx \
    --expectations    data/species_pathway_expectations.csv \
    --pathway-defs    data/pathway_definitions.csv \
    --output          runs/globdb_run/ko_grid_seq_subset.pdf \
    --title           "Sequenced MAGs – GlobDB-matched subset"
"""

import argparse
import re
import sys
from pathlib import Path

# ── import all plotting machinery from the existing script ────────────────────
sys.path.insert(0, str(Path(__file__).parent))
from plot_pathway_ko_grid import (
    load_pathway_defs, build_species_map, load_summary,
    build_grid, plot_grid, find_empty_pathways, load_ko_names,
    DEFAULT_IGNORE,
)

# GlobDB accession set (25 genomes in the globdb run)
GLOBDB_ACCESSIONS = {
    'GCA_023386795', 'GCF_000211415', 'GCF_000219215', 'GCF_000411455',
    'GCF_000697965', 'GCF_001541315', 'GCF_001591345', 'GCF_001898315',
    'GCF_002750975', 'GCF_002897035', 'GCF_003590645', 'GCF_003600685',
    'GCF_003729955', 'GCF_004330295', 'GCF_004514425', 'GCF_014217625',
    'GCF_014635885', 'GCF_021919345', 'GCF_024169245', 'GCF_024807945',
    'GCF_037023145', 'GCF_039727615', 'GCF_900116445', 'GCF_902506535',
    'GCF_943912325',
}


def load_globdb_genome_set(representatives_csv):
    """
    Parse representitives.csv and return the set of first-run genome IDs
    whose reference (stripping version suffix) is in GLOBDB_ACCESSIONS.
    """
    import csv
    kept = set()
    with open(representatives_csv, encoding='latin1') as f:
        for row in csv.DictReader(f):
            genome = row['User Genome'].replace('.fa_assembly', '').strip()
            ref = str(row.get('reference used', '') or '').strip()
            ref_base = ref.split('.')[0]   # strip version (e.g. .1, .2)
            if ref_base in GLOBDB_ACCESSIONS:
                kept.add(genome)
    return kept


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--representatives', required=True,
                    help='runs/globdb_run/representitives.csv')
    ap.add_argument('--summary',      required=True,
                    help='summary.xlsx from first run')
    ap.add_argument('--expectations', required=True,
                    help='data/species_pathway_expectations.csv')
    ap.add_argument('--pathway-defs', required=True,
                    help='data/pathway_definitions.csv')
    ap.add_argument('--names', action='store_true',
                    help='Show gene names instead of KO IDs')
    ap.add_argument('--genes-file', default='data/kegg_genes.csv')
    ap.add_argument('--ignore-pathway', action='append', default=[], metavar='PATHWAY')
    ap.add_argument('--include-all', action='store_true')
    ap.add_argument('--output',  default='ko_grid_globdb_subset.pdf')
    ap.add_argument('--title',   default='Sequenced MAGs – GlobDB-matched subset')
    args = ap.parse_args()

    ko_names = None
    if args.names:
        ko_names = load_ko_names(args.genes_file)

    print("Loading pathway definitions ...")
    pathway_defs = load_pathway_defs(args.pathway_defs)

    print("Loading expectations ...")
    species_map = build_species_map(args.expectations)
    all_expected = set()
    for pw_set in species_map.values():
        all_expected |= pw_set
    expected_pathway_order = [p for p in pathway_defs if p in all_expected]

    print("Determining GlobDB-matched genomes ...")
    keep = load_globdb_genome_set(args.representatives)
    print(f"  {len(keep)} genomes with a GlobDB reference: {sorted(keep)}")

    print("Loading summary Excel ...")
    genomes, species_labels, ko_hits, _ = load_summary(args.summary)
    print(f"  {len(genomes)} genomes in summary before filtering")

    genomes = [g for g in genomes if g in keep]
    print(f"  {len(genomes)} genomes after filtering")

    auto_ignore = [] if args.include_all else find_empty_pathways(genomes, ko_hits, expected_pathway_order)
    ignore = DEFAULT_IGNORE | set(args.ignore_pathway) | set(auto_ignore)
    if auto_ignore:
        print(f"  Auto-ignoring (absent in all): {', '.join(auto_ignore)}")
    expected_pathway_order = [p for p in expected_pathway_order if p not in ignore]
    print(f"  Pathways shown ({len(expected_pathway_order)}): {', '.join(expected_pathway_order)}")

    genomes = sorted(genomes, key=lambda g: re.sub(r'^s__', '', species_labels.get(g, g)).lower())

    print("Building grid ...")
    ko_columns, matrix, row_labels = build_grid(
        genomes, species_labels, ko_hits,
        pathway_defs, species_map, expected_pathway_order,
    )
    print(f"  Grid: {len(genomes)} genomes × {len(ko_columns)} KOs")

    print("Plotting ...")
    plot_grid(ko_columns, matrix, row_labels, args.title, args.output, ko_names=ko_names)


if __name__ == '__main__':
    main()
