"""
Small I/O helpers shared across phases.

Counts are stored as a genes x samples DataFrame; metadata as samples x
attributes. Everything round-trips through parquet in ``data/processed`` so
later phases never re-download or re-wrangle.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# Make the project root importable when scripts are run as ``python -m src.xxx``
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402


def save_counts(counts: pd.DataFrame, name: str = "counts.parquet") -> Path:
    """Persist a genes x samples count matrix."""
    path = config.processed_path(name)
    counts.to_parquet(path)
    print(f"[io] wrote counts {counts.shape} -> {path}")
    return path


def load_counts(name: str = "counts.parquet") -> pd.DataFrame:
    path = config.processed_path(name)
    return pd.read_parquet(path)


def save_metadata(meta: pd.DataFrame, name: str = "metadata.parquet") -> Path:
    """Persist a samples x attributes metadata table."""
    path = config.processed_path(name)
    meta.to_parquet(path)
    print(f"[io] wrote metadata {meta.shape} -> {path}")
    return path


def load_metadata(name: str = "metadata.parquet") -> pd.DataFrame:
    path = config.processed_path(name)
    return pd.read_parquet(path)


def save_table(df: pd.DataFrame, name: str, index: bool = True) -> Path:
    """Write a results table as CSV into results/tables."""
    config.TABLES.mkdir(parents=True, exist_ok=True)
    path = config.TABLES / f"{config.DATASET}_{name}"
    df.to_csv(path, index=index)
    print(f"[io] wrote table {df.shape} -> {path}")
    return path


def align_counts_metadata(
    counts: pd.DataFrame, meta: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Restrict both tables to the samples they share and put them in the same
    order. ``counts`` is genes x samples; ``meta`` is samples x attributes.
    """
    shared = [s for s in counts.columns if s in meta.index]
    dropped_c = set(counts.columns) - set(shared)
    dropped_m = set(meta.index) - set(shared)
    if dropped_c:
        print(f"[io] dropping {len(dropped_c)} count samples not in metadata")
    if dropped_m:
        print(f"[io] dropping {len(dropped_m)} metadata rows not in counts")
    counts = counts[shared]
    meta = meta.loc[shared]
    assert list(counts.columns) == list(meta.index), "sample alignment failed"
    print(f"[io] aligned to {len(shared)} shared samples")
    return counts, meta
