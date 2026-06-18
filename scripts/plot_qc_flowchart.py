"""QC flowchart: metagenome bins → GlobDB reference genomes."""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

# ── colour palette (colourblind-safe) ────────────────────────────────────────
C_MAIN   = '#2271B2'   # blue  – in analysis
C_FINAL  = '#359B73'   # green – final dataset
C_CONTAM = '#AA0000'   # red   – contamination removed
C_ASIDE1 = '#E69F00'   # amber – no GTDB assignment
C_ASIDE2 = '#7B5EA7'   # purple – not in GlobDB
C_TEXT   = 'white'
C_ARROW  = '#555555'

FIG_W, FIG_H = 17, 11.5

# ── column centres & box dimensions ──────────────────────────────────────────
X_MAIN = 4.2            # left spine
X_SIDE = 12.8           # right side-boxes
MW, MH = 5.8, 1.05     # main-box width / height
SW, MH_FINAL = 6.5, 1.3

# ── Y positions (main boxes determine row heights) ───────────────────────────
Y0 = 9.6    # 38 MAGs
Y1 = 7.0    # 34 MAGs
Y2 = 4.5    # 28 with GTDB
Y3 = 1.9    # 25 GlobDB (final)

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.axis('off')


# ─────────────────────────────────────────────────────────────────────────────
def rounded_box(ax, cx, cy, w, h, facecolor, zorder=3):
    p = FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                       boxstyle='round,pad=0.13',
                       facecolor=facecolor, edgecolor='white',
                       linewidth=2.2, zorder=zorder)
    ax.add_patch(p)


def main_box(cx, cy, label, sublabel=None, color=C_MAIN, h=MH):
    rounded_box(ax, cx, cy, MW, h, color)
    ty = cy + h * 0.13 if sublabel else cy
    ax.text(cx, ty, label, ha='center', va='center',
            fontsize=13, color=C_TEXT, fontweight='bold', zorder=4)
    if sublabel:
        ax.text(cx, cy - h * 0.26, sublabel, ha='center', va='center',
                fontsize=8, color=C_TEXT, zorder=4, style='italic')


def side_box(cx, cy, title, lines, color):
    n = len(lines)
    h = 0.55 + n * 0.36 + 0.15   # dynamic height
    rounded_box(ax, cx, cy, SW, h, color)
    top_y = cy + h / 2 - 0.32
    ax.text(cx, top_y, title, ha='center', va='center',
            fontsize=10.5, color=C_TEXT, fontweight='bold', zorder=4)
    for i, line in enumerate(lines):
        lx = cx - SW / 2 + 0.45
        ly = top_y - 0.38 - i * 0.35
        ax.text(lx, ly, line, ha='left', va='center',
                fontsize=8.3, color=C_TEXT, zorder=4, family='monospace')
    return h


def spine_arrow(x, y_from, y_to, label=None):
    ax.annotate('', xy=(x, y_to), xytext=(x, y_from),
                arrowprops=dict(arrowstyle='->', color=C_ARROW, lw=2.2), zorder=2)
    if label:
        ax.text(x - 0.35, (y_from + y_to) / 2, label,
                ha='right', va='center', fontsize=8.5,
                color='#444444', style='italic')


def branch_arrow(x_src, y, x_dst, color):
    ax.annotate('', xy=(x_dst, y), xytext=(x_src, y),
                arrowprops=dict(arrowstyle='->', color=color, lw=2.0), zorder=2)


# ─── Title ────────────────────────────────────────────────────────────────────
ax.text(FIG_W / 2, 11.0,
        'Metagenome bins → GlobDB reference genomes: QC funnel',
        ha='center', va='center', fontsize=15.5, fontweight='bold', color='#222222')

# ─── Main spine boxes ─────────────────────────────────────────────────────────
main_box(X_MAIN, Y0, '38 MAGs')
main_box(X_MAIN, Y1, '34 MAGs')
main_box(X_MAIN, Y2, '28 with GTDB assignment',
         sublabel='19 FastANI reference  ·  9 closest placement')
main_box(X_MAIN, Y3, '25 GlobDB reference genomes',
         sublabel='4 required taxonomy-based accession lookup',
         color=C_FINAL, h=MH_FINAL)

# ─── Spine arrows ─────────────────────────────────────────────────────────────
spine_arrow(X_MAIN, Y0 - MH / 2,  Y1 + MH / 2)
spine_arrow(X_MAIN, Y1 - MH / 2,  Y2 + MH / 2, label='GTDB-tk\nassignment')
spine_arrow(X_MAIN, Y2 - MH / 2,  Y3 + MH_FINAL / 2, label='GlobDB\nmatching')

# ─── Side box: contamination  (branches from 38-MAGs row) ────────────────────
h_c = side_box(
    X_SIDE, Y0,
    '4 removed  ·  contamination',
    ['SemiBin_32       Pseudomonas_E',
     'metabat2.21_sub  Methylobacterium',
     '425_sub          Shinella',
     '3854             Phyllobacterium'],
    C_CONTAM)
branch_arrow(X_MAIN + MW / 2, Y0, X_SIDE - SW / 2, C_CONTAM)

# ─── Side box: no GTDB assignment  (branches from 34-MAGs row) ───────────────
h_a1 = side_box(
    X_SIDE, Y1,
    '6 set aside  ·  no GTDB assignment',
    ['SemiBin_31   Pseudomonas_E   (genus only)',
     'SemiBin_27   unclassified',
     'SemiBin_17   Devosia_A       (genus only)',
     'metabat2.55  Variovorax      (genus only)',
     'metabat2.45  Scandinavium    (genus only)',
     '6176         Microbacterium  (genus only)'],
    C_ASIDE1)
branch_arrow(X_MAIN + MW / 2, Y1, X_SIDE - SW / 2, C_ASIDE1)

# ─── Side box: not in GlobDB  (branches from 28-with-GTDB row) ───────────────
h_a2 = side_box(
    X_SIDE, Y2,
    '4 set aside  ·  reference not in GlobDB  (lookup attempted)',
    ['SemiBin_11  Ochrobactrum',
     'SemiBin_23  Microbacterium',
     '6275        Agrobacterium',
     '803         Stenotrophomonas_A'],
    C_ASIDE2)
branch_arrow(X_MAIN + MW / 2, Y2, X_SIDE - SW / 2, C_ASIDE2)

# ─── Legend ───────────────────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(color=C_MAIN,   label='Genomes in analysis'),
    mpatches.Patch(color=C_FINAL,  label='Final dataset'),
    mpatches.Patch(color=C_CONTAM, label='Removed: contamination'),
    mpatches.Patch(color=C_ASIDE1, label='Set aside: no GTDB assignment'),
    mpatches.Patch(color=C_ASIDE2, label='Set aside: reference not in GlobDB'),
]
ax.legend(handles=legend_handles, loc='lower left',
          bbox_to_anchor=(0.005, 0.005), fontsize=9.5,
          framealpha=0.92, edgecolor='#cccccc')

plt.tight_layout(pad=0.4)
out_base = 'runs/globdb_run/qc_flowchart'
plt.savefig(f'{out_base}.pdf', bbox_inches='tight')
plt.savefig(f'{out_base}.png', bbox_inches='tight', dpi=180)
plt.close()
print('Saved', out_base + '.pdf / .png')
