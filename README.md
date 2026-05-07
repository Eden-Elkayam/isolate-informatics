# isolate-informatics

Snakemake pipeline for comparative genomic analysis of bacterial isolates. Takes Bakta-annotated genomes, runs KofamScan for KEGG ortholog annotation and Sourmash for GTDB taxonomy, then evaluates pathway presence using boolean gene logic and exports results as summary tables and heatmaps.

---

## Installation

**1. Create and activate the conda environment**

```bash
conda env create -f environment.yml
conda activate isolate_informatics
```

**2. Install the Python package**

```bash
pip install -e .
```

**3. Download required databases**

- **KofamScan** — profiles and KO list: https://www.genome.jp/ftp/db/kofam/
- **Sourmash GTDB** — database zip and lineages CSV: https://sourmash.readthedocs.io/en/latest/databases.html

By default the pipeline looks for these at `~/db/kofam` and `~/db/sourmash/`. You can override these paths in your run config.

---

## Running the pipeline

### 1. Create a folder for your run

```bash
mkdir -p runs/<run_name>
```

### 2. Copy and edit the config

```bash
cp config.yaml runs/<run_name>/config.yaml
```

Open `runs/<run_name>/config.yaml` and set `raw_dir` to the path containing your Bakta-annotated genomes. Each genome should be its own subdirectory holding the `.gb` (or `.gbff`) files.

```yaml
raw_dir: "/path/to/your/bakta_output"
```

Adjust any other settings as needed (KofamScan DB path, number of cores, output options).

### 3. Edit pathway expectations (optional)

`data/species_pathway_expectations.csv` lists known species–pathway associations used for validation. Add a row for any species in your run that you have prior knowledge about:

```
species,pathway,note,ref
My_organism,CBB,autotroph,Smith et al. 2020
```

The `species` field should match the GTDB-assigned taxonomy label for that genome. The `pathway` field must match a pathway name defined in `data/pathway_definitions.csv`.

### 4. Run

```bash
snakemake --configfile runs/<run_name>/config.yaml --cores 4
```

Results are written to the run folder. If `run_summary: true` is set in the config, a `summary.xlsx` file is produced. If `run_plot: true`, a KO presence heatmap is also generated.

---

## Key input files

| File | Description |
|---|---|
| `data/pathway_definitions.csv` | Boolean KO logic expressions defining each pathway |
| `data/kegg_genes.csv` | Human-readable gene names for KO IDs |
| `data/species_pathway_expectations.csv` | Known species–pathway associations for validation |
