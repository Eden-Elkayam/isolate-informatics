#!/usr/bin/env python3
"""
plot_comparison_ko_grid.py

Compare two summary Excel files and plot a KO × genome grid showing what changed.

Rows = union of genomes from both files (matched by genome ID).
        First-only genomes appear first, then genomes in both, then second-only.
Columns = KOs from pathway definitions (same layout as plot_pathway_ko_grid.py).

Visual encoding per cell:
  full color    — KO present in BOTH datasets for this genome
  hollow border — KO present in FIRST only (gene was lost / not called in second)
  sheer fill    — KO present in SECOND only (gene was gained / newly called)
  light grey    — KO absent in both

Usage
-----
python scripts/plot_comparison_ko_grid.py \\
    --first        runs/display/curated_meta.xlsx \\
    --second       runs/display/summary_glob.xlsx \\
    --pathway-defs data/pathway_definitions.csv \\
    --output       runs/display/ko_grid_comparison.pdf \\
    --title        "Comparison: curated MAGs vs GlobDB"
"""

import argparse
import re
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from plot_pathway_ko_grid import (
    DEFAULT_IGNORE,
    _PALETTE,
    _PATHWAY_COLORS,
    find_empty_pathways,
    load_ko_names,
    load_pathway_defs,
    load_summary,
)

# =============================================================================
# DIRECT RUN CONFIGURATION
# Fill these in and run the file directly from your IDE (no CLI needed).
# These are only used when no command-line arguments are passed.
# =============================================================================
DIRECT_RUN = dict(
    # ── inputs ────────────────────────────────────────────────────────────────
    first          = "runs/display/curated_meta.xlsx",   # baseline summary.xlsx
    second         = "runs/display/summary_glob.xlsx",   # comparison summary.xlsx
    pathway_defs   = "data/pathway_definitions.csv",

    # ── output ────────────────────────────────────────────────────────────────
    output         = "runs/display/ko_grid_comparison_withNames.pdf",  # .pdf or .png
    title          = "Comparison: curated MAGs vs GlobDB",

    # ── legend labels (how first/second are named in the legend) ──────────────
    first_label    = "curated MAGs",
    second_label   = "GlobDB",

    # ── pathway filtering ─────────────────────────────────────────────────────
    ignore_pathway = [          # list pathways to hide (exact names from pathway_defs)
        "H2 oxidation",
        "sulfide oxidation",
        "sulfite oxidation",
        "CO oxidation",
        "methanol oxidation",
        "nitrate reduction",
        "thiosulfate reduction",
    ],
    include_all    = False,     # True = show pathways absent in every genome too

    # ── column labels ─────────────────────────────────────────────────────────
    names          = True,     # True = show gene names instead of KO IDs
    genes_file     = "data/kegg_genes.csv",   # needed only when names=True
)
# =============================================================================

# Cell states
ABSENT         = 0
PRESENT_BOTH   = 1
PRESENT_FIRST  = 2
PRESENT_SECOND = 3

DARKEN_FACTOR = 0.6   # blend second-only cells with black to make them stronger
BORDER_LW     = 1.5


def _darken(color, factor=DARKEN_FACTOR):
    r, g, b = mcolors.to_rgb(color)
    return (r * factor, g * factor, b * factor)


def build_comparison_grid(
    genomes_first, genomes_second,
    species_first, species_second,
    ko_hits_first, ko_hits_second,
    pathway_defs, pathway_order,
):
    """
    Returns
    -------
    ko_columns   : list of (pathway, ko) in display order
    matrix       : ndarray (n_genomes × n_kos), values in {ABSENT, PRESENT_BOTH,
                   PRESENT_FIRST, PRESENT_SECOND}
    row_labels   : display label per row
    row_sources  : list of 'first', 'both', or 'second' per row
    """
    set_first  = set(genomes_first)
    set_second = set(genomes_second)

    first_only  = [g for g in genomes_first  if g not in set_second]
    in_both     = [g for g in genomes_first  if g in set_second]
    second_only = [g for g in genomes_second if g not in set_first]
    all_genomes = first_only + in_both + second_only

    # Merged species labels: first takes precedence, second fills gaps
    species_labels = {**species_second, **species_first}

    ko_columns = []
    for pathway in pathway_order:
        if pathway not in pathway_defs:
            continue
        for ko in pathway_defs[pathway]['kos']:
            ko_columns.append((pathway, ko))

    matrix      = np.zeros((len(all_genomes), len(ko_columns)), dtype=int)
    row_labels  = []
    row_sources = []

    for r, gid in enumerate(all_genomes):
        sp    = species_labels.get(gid, gid)
        label = re.sub(r'^s__', '', sp).strip()
        row_labels.append(label)

        in_f = gid in set_first
        in_s = gid in set_second
        row_sources.append('both' if (in_f and in_s) else ('first' if in_f else 'second'))

        for c, (pathway, ko) in enumerate(ko_columns):
            pf = ko in ko_hits_first.get((gid, pathway),  set())
            ps = ko in ko_hits_second.get((gid, pathway), set())
            if pf and ps:
                matrix[r, c] = PRESENT_BOTH
            elif pf:
                matrix[r, c] = PRESENT_FIRST
            elif ps:
                matrix[r, c] = PRESENT_SECOND
            # else ABSENT (default 0)

    return ko_columns, matrix, row_labels, row_sources, all_genomes


def plot_comparison_grid(
    ko_columns, matrix, row_labels, row_sources,
    title, output_path, ko_names=None,
    first_label='dataset 1', second_label='dataset 2',
):
    n_rows, n_cols = matrix.shape

    pathway_order, seen = [], set()
    for pw, _ in ko_columns:
        if pw not in seen:
            pathway_order.append(pw)
            seen.add(pw)

    _fallback = iter(c for c in _PALETTE if c not in _PATHWAY_COLORS.values())
    pw_color = {p: _PATHWAY_COLORS.get(p) or next(_fallback) for p in pathway_order}

    # ── geometry (same as plot_pathway_ko_grid) ───────────────────────────────
    cell_w     = 1.1
    cell_h     = 0.55
    label_w    = max((len(s) for s in row_labels), default=10) * 0.18 + 0.5
    ko_label_h = 5.0
    pw_bar_h   = 0.45
    pw_label_h = 3.0
    legend_h   = 1.0   # slightly taller to fit extra legend entries

    grid_w = n_cols * cell_w
    grid_h = n_rows * cell_h

    fig_w = max(14, label_w + grid_w + 0.5)
    fig_h = ko_label_h + grid_h + pw_bar_h + pw_label_h + legend_h + 0.3

    fig = plt.figure(figsize=(fig_w, fig_h))

    ax_x0 = label_w / fig_w
    ax_y0 = (ko_label_h + legend_h) / fig_h
    ax_w  = grid_w / fig_w
    ax_h  = grid_h / fig_h

    ax = fig.add_axes([ax_x0, ax_y0, ax_w, ax_h])
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.axis('off')

    # ── background stripe to distinguish source groups ────────────────────────
    group_start = 0
    current_src = row_sources[0] if row_sources else None
    stripe_colors = {'first': '#fffbe6', 'both': '#f0fff0', 'second': '#f0f4ff'}

    def _draw_stripe(src, start, end):
        row_y = n_rows - end
        h     = end - start
        ax.add_patch(plt.Rectangle(
            (-0.5, row_y - 0.5), n_cols + 0.5, h,
            facecolor=stripe_colors.get(src, 'white'),
            edgecolor='none', zorder=0,
        ))

    for r, src in enumerate(row_sources):
        if src != current_src:
            _draw_stripe(current_src, group_start, r)
            group_start  = r
            current_src  = src
    _draw_stripe(current_src, group_start, n_rows)

    # ── cells ─────────────────────────────────────────────────────────────────
    for r in range(n_rows):
        for c in range(n_cols):
            val   = matrix[r, c]
            pw    = ko_columns[c][0]
            row_y = n_rows - 1 - r
            color = pw_color[pw]

            if val == PRESENT_BOTH:
                # Shared — muted, not the interesting signal
                ax.add_patch(plt.Rectangle(
                    (c, row_y), 1, 1,
                    facecolor=color, edgecolor='white', linewidth=0.5, alpha=0.5,
                ))
            elif val == PRESENT_FIRST:
                # Lost in second — striped
                ax.add_patch(plt.Rectangle(
                    (c, row_y), 1, 1,
                    facecolor='white', edgecolor=color, linewidth=0.5,
                    hatch='////',
                ))
            elif val == PRESENT_SECOND:
                # New in second — full color, stands out
                ax.add_patch(plt.Rectangle(
                    (c, row_y), 1, 1,
                    facecolor=color, edgecolor='white', linewidth=0.5,
                ))
            else:
                ax.add_patch(plt.Rectangle(
                    (c, row_y), 1, 1,
                    facecolor='#f0f0f0', edgecolor='white', linewidth=0.5,
                ))

    # ── KO labels ─────────────────────────────────────────────────────────────
    for c, (_, ko) in enumerate(ko_columns):
        label = ko_names.get(ko, ko) if ko_names else ko
        fig.text(
            ax_x0 + (c + 0.5) * cell_w / fig_w,
            ax_y0 - 0.012,
            label,
            ha='right', va='top', fontsize=22, rotation=45,
            transform=fig.transFigure,
        )

    # ── row labels + source tag ────────────────────────────────────────────────
    for r, label in enumerate(row_labels):
        fig.text(
            ax_x0 + ax_w + 0.008,
            ax_y0 + (n_rows - 0.5 - r) * cell_h / fig_h,
            label,
            ha='left', va='center', fontsize=22,
            transform=fig.transFigure,
        )

    # ── pathway bars + labels ─────────────────────────────────────────────────
    pw_span = {}
    for c, (pw, _) in enumerate(ko_columns):
        if pw not in pw_span:
            pw_span[pw] = [c, c]
        pw_span[pw][1] = c

    bar_y0_fig = ax_y0 + ax_h + 0.008
    bar_h_fig  = pw_bar_h / fig_h
    lbl_y0_fig = bar_y0_fig + bar_h_fig + 0.005

    for pw in pathway_order:
        start, end = pw_span[pw]
        color = pw_color[pw]
        bx = ax_x0 + start * cell_w / fig_w
        bw = (end - start + 1) * cell_w / fig_w

        fig.add_artist(plt.Rectangle(
            (bx, bar_y0_fig), bw, bar_h_fig,
            facecolor=color, transform=fig.transFigure,
            figure=fig, clip_on=False,
        ))
        mid_x    = ax_x0 + (start + end + 1) / 2 * cell_w / fig_w
        pw_label = textwrap.fill(pw, width=12)
        fig.text(
            mid_x, lbl_y0_fig, pw_label,
            ha='center', va='bottom',
            fontsize=22, fontweight='bold', color=color,
            multialignment='center',
            transform=fig.transFigure,
        )
        if start > 0:
            sx = ax_x0 + start * cell_w / fig_w
            fig.add_artist(plt.Line2D(
                [sx, sx], [ax_y0, bar_y0_fig + bar_h_fig + 0.002],
                transform=fig.transFigure, figure=fig,
                color='#aaaaaa', linewidth=1.0, clip_on=False,
            ))

    # ── title ─────────────────────────────────────────────────────────────────
    title_y = bar_y0_fig + bar_h_fig + pw_label_h / fig_h + 0.01
    fig.text(
        ax_x0 + ax_w / 2, title_y,
        title, ha='center', va='bottom',
        fontsize=26, fontweight='bold',
        transform=fig.transFigure,
    )

    # ── legend ────────────────────────────────────────────────────────────────
    legend_handles = [
        mpatches.Patch(facecolor='#2ca02c', edgecolor='white', linewidth=0.5, alpha=0.5,
                       label=f'Present in both'),
        mpatches.Patch(facecolor='white',   edgecolor='#2ca02c', linewidth=0.5,
                       hatch='////', label=f'{first_label} only'),
        mpatches.Patch(facecolor='#2ca02c', edgecolor='white', linewidth=0.5,
                       label=f'{second_label} only (new)'),
        mpatches.Patch(facecolor='#f0f0f0', label='Absent in both'),
    ]
    fig.legend(
        handles=legend_handles,
        loc='lower left',
        bbox_to_anchor=(ax_x0, 0.01),
        bbox_transform=fig.transFigure,
        frameon=False, fontsize=20, ncol=2,
    )

    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    print(f"Saved: {output_path}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--first',        required=True, help='Baseline summary.xlsx')
    ap.add_argument('--second',       required=True, help='Comparison summary.xlsx')
    ap.add_argument('--pathway-defs', required=True, help='pathway_definitions.csv')
    ap.add_argument('--names',        action='store_true',
                    help='Show gene names instead of KO IDs where available')
    ap.add_argument('--genes-file',   default='data/kegg_genes.csv')
    ap.add_argument('--ignore-pathway', action='append', default=[], metavar='PATHWAY')
    ap.add_argument('--include-all',  action='store_true',
                    help='Include pathways absent in all genomes')
    ap.add_argument('--first-label',  default='dataset 1',
                    help='Display name for the first dataset (used in legend)')
    ap.add_argument('--second-label', default='dataset 2',
                    help='Display name for the second dataset (used in legend)')
    ap.add_argument('--output',  default='ko_grid_comparison.pdf')
    ap.add_argument('--title',   default='KO grid comparison')
    # If no CLI args, fall back to the DIRECT_RUN config block at the top
    if len(sys.argv) == 1:
        args = ap.parse_args([
            '--first',         DIRECT_RUN['first'],
            '--second',        DIRECT_RUN['second'],
            '--pathway-defs',  DIRECT_RUN['pathway_defs'],
            '--output',        DIRECT_RUN['output'],
            '--title',         DIRECT_RUN['title'],
            '--first-label',   DIRECT_RUN['first_label'],
            '--second-label',  DIRECT_RUN['second_label'],
            '--genes-file',    DIRECT_RUN['genes_file'],
            *(['--include-all'] if DIRECT_RUN['include_all'] else []),
            *(['--names']       if DIRECT_RUN['names']       else []),
            *[a for p in DIRECT_RUN['ignore_pathway'] for a in ('--ignore-pathway', p)],
        ])
    else:
        args = ap.parse_args()

    ko_names = None
    if args.names:
        print(f"Loading KO names from {args.genes_file} ...")
        ko_names = load_ko_names(args.genes_file)

    print("Loading pathway definitions ...")
    pathway_defs = load_pathway_defs(args.pathway_defs)
    all_pathways = list(pathway_defs.keys())

    print("Loading first summary ...")
    genomes_first, species_first, ko_hits_first, pathways_in_first = load_summary(args.first)
    print(f"  {len(genomes_first)} genomes")

    print("Loading second summary ...")
    genomes_second, species_second, ko_hits_second, pathways_in_second = load_summary(args.second)
    print(f"  {len(genomes_second)} genomes")

    # Pathway order: union of what appears in either file, filtered to definitions
    pathway_union = dict.fromkeys(
        p for p in all_pathways
        if p in pathways_in_first or p in pathways_in_second
    )
    pathway_order = list(pathway_union)

    # Auto-ignore pathways absent in all genomes across both files
    all_genomes_combined = list(dict.fromkeys(genomes_first + genomes_second))
    combined_hits = {**ko_hits_first, **ko_hits_second}
    if not args.include_all:
        auto_ignore = find_empty_pathways(all_genomes_combined, combined_hits, pathway_order)
        if auto_ignore:
            print(f"  Auto-ignoring (absent in all genomes): {', '.join(auto_ignore)}")
    else:
        auto_ignore = []

    ignore = DEFAULT_IGNORE | set(args.ignore_pathway) | set(auto_ignore)
    pathway_order = [p for p in pathway_order if p not in ignore]
    print(f"  Pathways shown ({len(pathway_order)}): {', '.join(pathway_order)}")

    overlap = set(genomes_first) & set(genomes_second)
    print(f"  Genomes: {len(genomes_first)} first, {len(genomes_second)} second, "
          f"{len(overlap)} in both")

    print("Building comparison grid ...")
    ko_columns, matrix, row_labels, row_sources, all_genomes = build_comparison_grid(
        genomes_first, genomes_second,
        species_first, species_second,
        ko_hits_first, ko_hits_second,
        pathway_defs, pathway_order,
    )
    print(f"  Grid: {len(all_genomes)} genomes × {len(ko_columns)} KOs")

    print("Plotting ...")
    plot_comparison_grid(
        ko_columns, matrix, row_labels, row_sources,
        args.title, args.output, ko_names=ko_names,
        first_label=args.first_label, second_label=args.second_label,
    )


if __name__ == '__main__':
    main()
