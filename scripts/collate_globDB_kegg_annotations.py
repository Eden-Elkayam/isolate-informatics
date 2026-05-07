#!/usr/bin/env python3
"""
Collate KEGG annotations from globDB protein annotation files.

This script walks an input directory, opens gzipped tab-delimited annotation files
(.tsv.gz), extracts KEGG annotation fields, and counts KEGG IDs
per genome (top-level subdirectory name).

Usage:
  python scripts/collate_globDB_kegg_annotations.py \
    --input ~/Downloads/globdb_r226_protein_annotations \
    --output output/globdb_kegg_counts.tsv
"""

from __future__ import annotations

import argparse
import os
import threading
from concurrent.futures import ProcessPoolExecutor, as_completed
from queue import Queue
from tqdm import tqdm
from collections import Counter, defaultdict
from typing import Dict, Iterable

import pandas as pd


def load_valid_genomes(metadata_path: str, min_completeness: float | None = None, max_contamination: float | None = None) -> set[str]:
    """Load set of valid genome IDs from metadata file based on quality filters.
    
    Args:
        metadata_path: Path to metadata TSV file
        min_completeness: Minimum checkm_completeness threshold (0-100)
        max_contamination: Maximum checkm_contamination threshold (0-100)
    
    Returns:
        Set of valid genome IDs
    """
    df = pd.read_csv(metadata_path, sep="\t")
    total_genomes = len(df)
    
    # Check required columns exist
    completeness_col = "checkm2_completeness"
    contamination_col = "checkm2_contamination"
    if min_completeness is not None and completeness_col not in df.columns:
        raise ValueError(f"{completeness_col} column not found in {metadata_path}")
    if max_contamination is not None and contamination_col not in df.columns:
        raise ValueError(f"{contamination_col} column not found in {metadata_path}")
    
    # Apply filters
    mask = pd.Series([True] * len(df))
    if min_completeness is not None:
        mask &= df[completeness_col] >= min_completeness
    if max_contamination is not None:
        mask &= df[contamination_col] <= max_contamination
    
    valid_genomes = set(df.loc[mask].iloc[:, 0].astype(str))  # First column assumed to be genome ID
    filtered_out = total_genomes - len(valid_genomes)
    fraction_filtered = filtered_out / total_genomes if total_genomes > 0 else 0
    
    print(f"Quality filter results from {metadata_path}:")
    print(f"  Total genomes: {total_genomes}")
    print(f"  Passed filter: {len(valid_genomes)} ({100 * (1 - fraction_filtered):.1f}%)")
    print(f"  Filtered out: {filtered_out} ({100 * fraction_filtered:.1f}%)")
    print(f"  Filters: {completeness_col}>={min_completeness}, {contamination_col}<={max_contamination}")
    
    return valid_genomes


def iter_annotation_files(input_dir: str) -> Iterable[str]:
    """Recursively yield paths to all .tsv.gz annotation files in input_dir."""
    for root, _, files in os.walk(input_dir, followlinks=True):
        if root != input_dir:
            print(f"Entering subdirectory {root}")
        for fname in files:
            if fname.startswith("."):
                continue
            if fname.endswith(".tsv.gz"):
                yield os.path.join(root, fname)


def genome_from_path(input_dir: str, archive_path: str) -> str:
    """Extract genome identifier from filename by stripping .tsv.gz and _annotations suffix."""
    base = os.path.basename(archive_path)
    if base.endswith(".tsv.gz"):
        base = base[: -len(".tsv.gz")]
    if base.endswith("_annotations"):
        base = base[: -len("_annotations")]
    return base or "root"


def count_kegg_in_file(path: str) -> Dict[str, int]:
    """Read gzipped TSV, filter for KOfam source, and return KEGG ID counts."""
    df = pd.read_csv(path, sep="\t", compression="gzip", dtype=str)
    if df.empty:
        return {}

    mask = df["source"].astype(str).str.strip().eq("KOfam")
    accessions = df.loc[mask, "accession"].astype(str).str.strip()
    return accessions.value_counts().to_dict()


def _process_single_file(archive_path: str, input_dir: str) -> tuple[str, Dict[str, int]]:
    """Helper for parallel processing: returns (genome, kegg_counts) tuple."""
    genome = genome_from_path(input_dir, archive_path)
    file_counts = count_kegg_in_file(archive_path)
    return genome, file_counts



def _batch_writer_thread(output_path: str, queue: Queue, batch_size: int = 10000, min_itemsize: Dict[str, int] | None = None) -> None:
    """Background thread that writes batches of rows to HDF5.
    
    Receives tuples of (genome, kegg_id, count) from queue and writes in batches.
    Stops when queue receives None sentinel.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True) if os.path.dirname(output_path) else None
    rows = []
    first_write = True
    
    while True:
        item = queue.get()
        if item is None:  # Sentinel: stop signal
            # Write final batch
            if rows:
                df = pd.DataFrame(rows)
                df.to_hdf(output_path, key="kegg_counts", mode="a", append=True, format="table", min_itemsize=min_itemsize)
            break
        
        rows.append(item)
        
        # Write batch if full
        if len(rows) >= batch_size:
            df = pd.DataFrame(rows)
            mode = "w" if first_write else "a"
            df.to_hdf(output_path, key="kegg_counts", mode=mode, append=not first_write, format="table", min_itemsize=min_itemsize)
            first_write = False
            rows = []


def write_counts_hdf5_streaming(counts: Dict[str, Counter], output_path: str, batch_size: int = 10000) -> None:
    """Write counts to HDF5 by streaming rows through a queue.
    
    Uses a background writer thread to batch writes, reducing memory overhead.
    """
    print(f"Writing output to {output_path}")
    
    # Pre-compute max string lengths to avoid HDF5 column size validation errors
    max_genome_len = max(len(g) for g in counts.keys()) if counts else 1
    max_kegg_len = max(max(len(k) for k in genome_counts.keys()) for genome_counts in counts.values()) if counts else 1
    
    queue: Queue = Queue(maxsize=batch_size * 2)
    
    # Start writer thread with min_itemsize dict
    min_itemsize = {"genome": max_genome_len, "kegg_id": max_kegg_len}
    writer = threading.Thread(
        target=_batch_writer_thread,
        args=(output_path, queue, batch_size, min_itemsize),
        daemon=False
    )
    writer.start()
    
    # Submit rows from dict to queue
    total_rows = 0
    for genome in sorted(counts.keys()):
        for kegg_id, count in counts[genome].most_common():
            queue.put({"genome": genome, "kegg_id": kegg_id, "count": count})
            total_rows += 1
    
    # Signal end of data
    queue.put(None)
    writer.join()
    
    print(f"Wrote {total_rows} rows to {output_path}")


def collate_kegg_counts(input_dir: str, num_workers: int | None = None, valid_genomes: set[str] | None = None) -> Dict[str, Counter]:
    """Process all annotation files in parallel and aggregate KEGG counts by genome.
    
    Args:
        input_dir: Directory containing annotation files
        num_workers: Number of parallel workers
        valid_genomes: Optional set of genome IDs to include (others will be skipped)
    """
    counts: Dict[str, Counter] = defaultdict(Counter)

    archive_paths = list(iter_annotation_files(input_dir))
    # TEMPORARY: Limit to 100 files for testing
    #archive_paths = archive_paths[:100]

    print(f"Found {len(archive_paths)} annotation files to process.")
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(_process_single_file, path, input_dir): path for path in archive_paths}
        processed = 0
        skipped = 0
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing archives"):
            genome, file_counts = future.result()
            
            # Filter by valid genomes if provided
            if valid_genomes is not None and genome not in valid_genomes:
                skipped += 1
                continue
            
            processed += 1
            for kegg_id, count in file_counts.items():
                counts[genome][kegg_id] += count
        
        if valid_genomes is not None:
            print(f"Processed {processed} genomes, skipped {skipped} (filtered by quality thresholds)")

    return counts


def build_arg_parser() -> argparse.ArgumentParser:
    """Build and return argument parser for command-line interface."""
    parser = argparse.ArgumentParser(description="Collate KEGG annotation counts per genome from globDB annotations.")
    parser.add_argument("-i", "--input", required=True, help="Path to globDB protein annotations directory.")
    parser.add_argument("-o", "--output", required=True, help="Output HDF5 path.")
    parser.add_argument("-n", "--num-workers", type=int, default=None, help="Number of parallel workers (default: number of CPU cores).")
    parser.add_argument("-m", "--metadata", default=None, help="Path to metadata TSV file with genome quality metrics.")
    parser.add_argument("--min-completeness", type=float, default=None, help="Minimum checkm_completeness threshold (0-100).")
    parser.add_argument("--max-contamination", type=float, default=None, help="Maximum checkm_contamination threshold (0-100).")
    return parser


def main() -> None:
    """Parse arguments, collate KEGG counts, and write output."""
    parser = build_arg_parser()
    args = parser.parse_args()
    
    # Validate metadata arguments
    if (args.min_completeness is not None or args.max_contamination is not None) and args.metadata is None:
        parser.error("--metadata is required when --min-completeness or --max-contamination is specified")
    
    # Load valid genomes if filtering requested
    valid_genomes = None
    if args.metadata is not None:
        valid_genomes = load_valid_genomes(args.metadata, args.min_completeness, args.max_contamination)

    counts = collate_kegg_counts(args.input, num_workers=args.num_workers, valid_genomes=valid_genomes)
    write_counts_hdf5_streaming(counts, args.output)


if __name__ == "__main__":
    main()
