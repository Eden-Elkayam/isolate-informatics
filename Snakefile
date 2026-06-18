"""
Isolate genome annotation pipeline.

Preprocesses Bakta-annotated isolate genomes through KofamScan + Sourmash
taxonomy, then evaluates pathway presence using boolean logic.

Usage:
  conda activate c1_informatics
  snakemake --configfile runs/<run_name>/config.yaml --cores 4
"""

import sys
from pathlib import Path

PYTHON = sys.executable

configfile: "config.yaml" if Path("config.yaml").exists() else {}

# ── Config defaults ────────────────────────────────────────────────────────────
config.setdefault("raw_dir", "")
config.setdefault("gb_ext", "gb")
config.setdefault("kofam_db", "~/db/kofam")
config.setdefault("kofamscan_cpus", 4)
config.setdefault("sourmash_db", "~/db/sourmash/gtdb-rs214-reps.k31.zip")
config.setdefault("sourmash_lineages", "~/db/sourmash/taxonomy/gtdb-rs214.lineages.csv")
# Auto-derive output_dir from the run configfile's directory so outputs are
# always co-located with the config that produced them.
if not config.get("output_dir") or config["output_dir"] in ("output", "null", None, ""):
    _run_cf = next(
        (cf for cf in workflow.configfiles
         if Path(cf).resolve() != Path("config.yaml").resolve()),
        None,
    )
    config["output_dir"] = str(Path(_run_cf).parent) if _run_cf else "output"
config.setdefault("num_workers", 4)
config.setdefault("min_annotations_per_genome", 10)
config.setdefault("min_genomes_per_annotation", 1)
config.setdefault("kegg_genes_file", None)
config.setdefault("pathway_definitions_file", None)
config.setdefault("run_summary", False)
config.setdefault("run_plot", False)
config.setdefault("plot_title", "Pathway KO Presence Grid")
config.setdefault("plot_format", "pdf")
config.setdefault("plot_ignore_pathways", "")
config.setdefault("plot_include_empty", False)
config.setdefault("plot_show_names", False)

# ── Input mode ─────────────────────────────────────────────────────────────────
FNA_MODE = config["gb_ext"] == "fna"
FAA_MODE = config["gb_ext"] == "faa"
# Taxonomy requires nucleotide sequences; available in Bakta and FNA modes.
config.setdefault("run_taxonomy", not FAA_MODE)

# ── Paths ──────────────────────────────────────────────────────────────────────
RAW_DIR = Path(config["raw_dir"]).expanduser()
OUTPUT_DIR = Path(config["output_dir"])
KEGG_GENES_FILE = Path(config["kegg_genes_file"]) if config["kegg_genes_file"] else None
PATHWAY_DEFINITIONS_FILE = Path(config["pathway_definitions_file"]) if config["pathway_definitions_file"] else None

KEGG_DIR = OUTPUT_DIR / "kegg"
MERGED_DIR = KEGG_DIR / "merged"
BAKTA_INPUT_DIR = KEGG_DIR / "bakta_input"
KO_DIR = KEGG_DIR / "kofamscan"
KOFAMSCAN_INPUT_DIR = KEGG_DIR / "kofamscan_input"

PHYLO_DIR = OUTPUT_DIR / "phylogeny"
FNA_DIR = PHYLO_DIR / "fna"
TAXONOMY_DIR = PHYLO_DIR / "taxonomy"
TAXONOMY_SUMMARY = TAXONOMY_DIR / "taxonomy_summary.csv"

PATHWAY_DIR = OUTPUT_DIR / "pathways"
COUNTS_FILE = PATHWAY_DIR / "kegg_counts.h5"
MATRIX_FILE = PATHWAY_DIR / "kegg_by_genome_matrix.h5"
FILTERED_MATRIX_FILE = PATHWAY_DIR / "kegg_subset_matrix.h5"
PATHWAY_VALIDATION_FILE = PATHWAY_DIR / "pathway_gene_references_validated.txt"
PATHWAY_MATRIX_FILE = PATHWAY_DIR / "pathway_presence_matrix.h5"
SUMMARY_XLSX = OUTPUT_DIR / "summary.xlsx"
KO_GRID_FILE = OUTPUT_DIR / f"ko_grid.{config['plot_format']}"

if FNA_MODE:
    # raw_dir is a flat directory of assembled .fna genomes.
    # Prodigal predicts proteins; .fna files are used directly for taxonomy.
    FAA_DIR = KEGG_DIR / "faa"
    GENOMES = sorted(f.stem for f in RAW_DIR.glob("*.fna")) if RAW_DIR.exists() else []
    C1_INPUT_DIR = KOFAMSCAN_INPUT_DIR
elif FAA_MODE:
    # raw_dir is a flat directory of pre-computed prodigal .faa files.
    # Genome IDs are file stems (e.g. "2445.fa" from "2445.fa.faa").
    FAA_DIR = RAW_DIR
    GENOMES = sorted(f.stem for f in RAW_DIR.glob("*.faa")) if RAW_DIR.exists() else []
    C1_INPUT_DIR = KOFAMSCAN_INPUT_DIR
else:
    # Bakta mode: FAA files are extracted from merged .gb files.
    FAA_DIR = KEGG_DIR / "faa"
    GENOMES = sorted(f"{d.name}.merged" for d in RAW_DIR.iterdir() if d.is_dir()) if RAW_DIR.exists() else []
    C1_INPUT_DIR = KEGG_DIR / "c1_input"


# ── Validation ─────────────────────────────────────────────────────────────────
def validate_config():
    errors = []
    if not config.get("raw_dir"):
        errors.append("raw_dir is required")
    elif not RAW_DIR.exists():
        errors.append(f"raw_dir does not exist: {RAW_DIR}")
    if KEGG_GENES_FILE and not KEGG_GENES_FILE.exists():
        errors.append(f"kegg_genes_file does not exist: {KEGG_GENES_FILE}")
    if PATHWAY_DEFINITIONS_FILE and not PATHWAY_DEFINITIONS_FILE.exists():
        errors.append(f"pathway_definitions_file does not exist: {PATHWAY_DEFINITIONS_FILE}")
    if errors:
        raise ValueError("Configuration validation failed:\n  " + "\n  ".join(errors))

validate_config()


# ── Helpers ────────────────────────────────────────────────────────────────────
def build_ignore_args():
    ignore_str = config.get("plot_ignore_pathways", "")
    if not ignore_str:
        return ""
    pathways = [p.strip() for p in str(ignore_str).split(",") if p.strip()]
    return " ".join(f'--ignore-pathway "{p}"' for p in pathways)


# ── Rule all ───────────────────────────────────────────────────────────────────
rule all:
    input:
        str(C1_INPUT_DIR),
        str(TAXONOMY_SUMMARY) if config.get("run_taxonomy") else [],
        MATRIX_FILE,
        PATHWAY_DIR / "annotations_per_genome_histogram.png",
        PATHWAY_DIR / "genomes_per_kegg_histogram.png",
        FILTERED_MATRIX_FILE if KEGG_GENES_FILE else [],
        PATHWAY_MATRIX_FILE if PATHWAY_DEFINITIONS_FILE else [],
        SUMMARY_XLSX if config.get("run_summary") else [],
        str(KO_GRID_FILE) if config.get("run_plot") else [],


# ── KEGG preprocessing ─────────────────────────────────────────────────────────

rule merge_gbk:
    """Merge per-contig .gb files into one file per genome."""
    input:
        str(RAW_DIR),
    output:
        directory(MERGED_DIR),
    shell:
        "{PYTHON} scripts/merge_gbk.py {input} --out {output} --ext " + str(config["gb_ext"])


rule extract_proteins:
    """Convert merged .gb files to .faa protein FASTAs for KofamScan."""
    input:
        MERGED_DIR,
    output:
        directory(FAA_DIR),
    shell:
        "{PYTHON} scripts/gb_to_faa.py {input} --out {output} --ext " + str(config["gb_ext"])


rule extract_bakta_kegg:
    """Extract Bakta KEGG annotations from merged .gb files."""
    input:
        MERGED_DIR,
    output:
        directory(BAKTA_INPUT_DIR),
    shell:
        "{PYTHON} scripts/bakta_to_pipeline_input.py --input-dir {input} --output-dir {output}"


if FNA_MODE:
    rule prodigal:
        """Predict protein-coding genes from an assembled genome using Prodigal."""
        input:
            fna=str(RAW_DIR / "{genome}.fna"),
        output:
            faa=str(FAA_DIR / "{genome}.faa"),
        shell:
            "prodigal -i {input.fna} -a {output.faa} -p single -f gff -o /dev/null"


rule kofamscan:
    """Assign KEGG KOs to one genome via KofamScan."""
    input:
        faa=str(FAA_DIR / "{genome}.faa"),
    output:
        tsv=str(KO_DIR / "{genome}.ko.tsv"),
    params:
        profiles=str(Path(config["kofam_db"]).expanduser() / "profiles"),
        ko_list=str(Path(config["kofam_db"]).expanduser() / "ko_list"),
        cpus=config["kofamscan_cpus"],
        tmp_dir=str(KO_DIR / "{genome}_tmp"),
    shell:
        """
        exec_annotation \
            -p {params.profiles} \
            -k {params.ko_list} \
            --cpu {params.cpus} \
            -f detail-tsv \
            --tmp-dir {params.tmp_dir} \
            -o {output.tsv} \
            {input.faa}
        rm -rf {params.tmp_dir}
        """


rule convert_kofamscan:
    """Convert KofamScan outputs to c1_informatics input format."""
    input:
        expand(str(KO_DIR / "{genome}.ko.tsv"), genome=GENOMES),
    output:
        directory(KOFAMSCAN_INPUT_DIR),
    shell:
        "{PYTHON} scripts/kofamscan_to_pipeline_input.py --input-dir " + str(KO_DIR) + " --output-dir {output}"


rule merge_kegg_annotations:
    """Merge Bakta and KofamScan KEGG annotations (union, deduplicated)."""
    input:
        bakta=BAKTA_INPUT_DIR,
        kofamscan=KOFAMSCAN_INPUT_DIR,
    output:
        directory(C1_INPUT_DIR),
    shell:
        """
        {PYTHON} scripts/merge_kegg_inputs.py \
            --bakta-dir {input.bakta} \
            --kofamscan-dir {input.kofamscan} \
            --output-dir {output}
        """


rule audit_kegg:
    """Audit KEGG annotation coverage for target pathways."""
    input:
        MERGED_DIR,
    output:
        directory(OUTPUT_DIR / "kegg_audit"),
    shell:
        "{PYTHON} scripts/find_missing_kegg.py {input} --out {output} --ext " + str(config["gb_ext"])


# ── Phylogeny ──────────────────────────────────────────────────────────────────

rule extract_nucleotides:
    """Extract nucleotide FASTAs from merged .gb files for Sourmash."""
    input:
        MERGED_DIR,
    output:
        directory(FNA_DIR),
    shell:
        "{PYTHON} scripts/gb_to_fna.py {input} --out {output} --ext " + str(config["gb_ext"])


rule sourmash_taxonomy:
    """Identify GTDB species for each genome via Sourmash + GTDB database."""
    input:
        # FNA mode: assembled .fna files are in raw_dir; no extraction needed.
        RAW_DIR if FNA_MODE else FNA_DIR,
    output:
        lineage=str(TAXONOMY_DIR / "gtdb_taxonomy.lineage.csv"),
        summary=str(TAXONOMY_SUMMARY),
    params:
        db=str(Path(config["sourmash_db"]).expanduser()),
        lineages=str(Path(config["sourmash_lineages"]).expanduser()),
        output_dir=str(TAXONOMY_DIR),
    shell:
        """
        {PYTHON} scripts/sourmash_taxonomy.py \
            --fna-dir {input} \
            --db {params.db} \
            --lineages {params.lineages} \
            --output-dir {params.output_dir}
        """


# ── Pathway analysis ───────────────────────────────────────────────────────────

rule collate_kegg_counts:
    """Aggregate KEGG annotations into counts per genome."""
    input:
        str(C1_INPUT_DIR),
    output:
        counts=COUNTS_FILE,
    params:
        num_workers=config["num_workers"],
    shell:
        """
        {PYTHON} scripts/collate_globDB_kegg_annotations.py \
            -i {input} \
            -o {output.counts} \
            -n {params.num_workers}
        """


rule create_kegg_matrix:
    """Filter KEGG counts and create genome x annotation matrix."""
    input:
        counts=COUNTS_FILE,
    output:
        matrix=MATRIX_FILE,
        hist1=PATHWAY_DIR / "annotations_per_genome_histogram.png",
        hist2=PATHWAY_DIR / "genomes_per_kegg_histogram.png",
    params:
        min_annotations=config["min_annotations_per_genome"],
        min_genomes=config["min_genomes_per_annotation"],
    shell:
        """
        {PYTHON} scripts/kegg_counts_to_matrix.py \
          -i {input.counts} \
          -o {output.matrix} \
          -a {params.min_annotations} \
          -g {params.min_genomes}
        """


rule filter_kegg_matrix:
    """Filter matrix to subset of target KEGG genes."""
    input:
        matrix=MATRIX_FILE,
        genes=KEGG_GENES_FILE,
    output:
        filtered=FILTERED_MATRIX_FILE,
    shell:
        """
        {PYTHON} scripts/filter_kegg_matrix.py \
          -i {input.matrix} \
          -g {input.genes} \
          -o {output.filtered}
        """


rule validate_pathway_gene_references:
    """Ensure pathway boolean expressions only reference known KEGG genes."""
    input:
        pathways=PATHWAY_DEFINITIONS_FILE,
        genes=KEGG_GENES_FILE,
    output:
        validated=PATHWAY_VALIDATION_FILE,
    shell:
        """
        {PYTHON} scripts/validate_pathway_gene_references.py \
            --pathways {input.pathways} \
            --genes {input.genes}
        touch {output.validated}
        """


rule evaluate_pathway_logic:
    """Evaluate pathway presence using boolean logic over KEGG genes."""
    input:
        matrix=FILTERED_MATRIX_FILE,
        pathways=PATHWAY_DEFINITIONS_FILE,
        genes=KEGG_GENES_FILE,
        validated=PATHWAY_VALIDATION_FILE,
    output:
        pathway_matrix=PATHWAY_MATRIX_FILE,
    shell:
        """
        {PYTHON} scripts/evaluate_pathway_logic.py \
          -m {input.matrix} \
          -p {input.pathways} \
          -g {input.genes} \
          -o {output.pathway_matrix}
        """


# ── Summary and plot ───────────────────────────────────────────────────────────

rule export_excel:
    """Export pathway and gene summary tables to Excel."""
    input:
        pathway_matrix=PATHWAY_MATRIX_FILE,
        gene_matrix=FILTERED_MATRIX_FILE if KEGG_GENES_FILE else [],
        taxonomy=[str(TAXONOMY_SUMMARY)] if config.get("run_taxonomy") else [],
    output:
        xlsx=SUMMARY_XLSX,
    params:
        output_dir=str(PATHWAY_DIR),
        taxonomy_arg=f"--taxonomy {TAXONOMY_SUMMARY}" if config.get("run_taxonomy") else "",
    shell:
        """
        {PYTHON} scripts/export_excel.py \
            --output-dir {params.output_dir} \
            {params.taxonomy_arg} \
            --output {output.xlsx}
        """


rule plot_ko_grid:
    """Plot KO presence grid from summary Excel."""
    input:
        summary=SUMMARY_XLSX,
        expectations="data/species_pathway_expectations.csv",
        pathway_defs=PATHWAY_DEFINITIONS_FILE,
    output:
        plot=str(KO_GRID_FILE),
    params:
        title=config["plot_title"],
        ignore_args=build_ignore_args(),
        include_empty="--include-all" if config.get("plot_include_empty") else "",
        show_names="--names --genes-file data/kegg_genes.csv" if config.get("plot_show_names") else "",
    shell:
        """
        {PYTHON} scripts/plot_pathway_ko_grid.py \
            --expectations {input.expectations} \
            --summary {input.summary} \
            --pathway-defs {input.pathway_defs} \
            --output {output.plot} \
            --title "{params.title}" \
            {params.ignore_args} \
            {params.include_empty} \
            {params.show_names}
        """


rule clean:
    """Remove all generated outputs."""
    run:
        import shutil
        shutil.rmtree(OUTPUT_DIR, ignore_errors=True)
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        print(f"Removed {OUTPUT_DIR}")
