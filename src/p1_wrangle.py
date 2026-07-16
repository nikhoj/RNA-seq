"""
Phase 1b - Data wrangling + first QC visuals.

Takes the raw counts/metadata cached by p1_download and produces analysis-ready
tables:
  - subset/relabel samples per the dataset DESIGN spec
  - parse condition (and time) into tidy ordered factors
  - align sample IDs between counts and metadata
  - filter low-count genes
  - QC visuals: library-size bars, detected-genes-per-sample, count densities

Run:  RNASEQ_DATASET=airway python -m src.p1_wrangle
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402
from src import io_utils, theme  # noqa: E402,F401


# --------------------------------------------------------------------------- #
# Factor parsing
# --------------------------------------------------------------------------- #
def _tidy_condition(meta: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """Build a clean, ordered ``condition`` column with baseline first."""
    col = spec["condition_col"]
    if col not in meta.columns:
        raise KeyError(
            f"condition_col '{col}' not in metadata columns {list(meta.columns)}"
        )
    meta = meta.copy()
    if spec.get("condition_keep"):
        meta = meta[meta[col].isin(spec["condition_keep"])]

    # normalise to short lowercase labels
    meta["condition"] = (
        meta[col].astype(str).str.strip().str.lower().str.replace(" ", "_")
    )
    baseline = str(spec["condition_baseline"]).lower()
    levels = [baseline] + sorted(l for l in meta["condition"].unique() if l != baseline)
    meta["condition"] = pd.Categorical(meta["condition"], categories=levels, ordered=True)
    return meta


def _tidy_time(meta: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """Parse a numeric ``time`` column if the dataset has a time axis."""
    tcol = spec.get("time_col")
    if not tcol or tcol not in meta.columns:
        return meta
    meta = meta.copy()
    # extract the first number found (handles "day 4", "4d", "T4", etc.)
    nums = meta[tcol].astype(str).str.extract(r"(-?\d+\.?\d*)")[0].astype(float)
    meta["time"] = nums
    meta["time_factor"] = pd.Categorical(
        meta["time"].astype("Int64").astype(str),
        categories=[str(int(t)) for t in sorted(meta["time"].dropna().unique())],
        ordered=True,
    )
    return meta


# --------------------------------------------------------------------------- #
# Gene filtering
# --------------------------------------------------------------------------- #
def filter_low_counts(counts: pd.DataFrame) -> pd.DataFrame:
    """Keep genes with >= MIN_COUNT in >= MIN_SAMPLES samples."""
    keep = (counts >= config.MIN_COUNT).sum(axis=1) >= config.MIN_SAMPLES
    filtered = counts[keep]
    print(
        f"[p1] gene filter: {keep.sum():,}/{len(counts):,} genes kept "
        f"(>= {config.MIN_COUNT} counts in >= {config.MIN_SAMPLES} samples)"
    )
    return filtered


# --------------------------------------------------------------------------- #
# QC visuals
# --------------------------------------------------------------------------- #
def qc_plots(counts: pd.DataFrame, meta: pd.DataFrame) -> None:
    from plotnine import (
        aes, geom_bar, geom_col, geom_density, ggplot, labs,
        scale_color_manual, scale_fill_manual, scale_x_log10,
    )

    outdir = config.fig_dir("phase1_qc")
    pal = theme.OKABE_ITO

    # --- library size per sample -------------------------------------------
    lib = pd.DataFrame({
        "sample": counts.columns,
        "library_size": counts.sum(axis=0).values / 1e6,
        "condition": meta["condition"].astype(str).values,
    })
    p = (
        ggplot(lib, aes("sample", "library_size", fill="condition"))
        + geom_col()
        + scale_fill_manual(values=pal)
        + labs(title="Library size per sample", x="", y="Total counts (millions)")
        + theme.GGPLOT_THEME
    )
    _rotate_x(p, counts.shape[1]).save(outdir / "library_size.png", verbose=False)

    # --- detected genes per sample -----------------------------------------
    det = pd.DataFrame({
        "sample": counts.columns,
        "detected_genes": (counts > 0).sum(axis=0).values,
        "condition": meta["condition"].astype(str).values,
    })
    p = (
        ggplot(det, aes("sample", "detected_genes", fill="condition"))
        + geom_col()
        + scale_fill_manual(values=pal)
        + labs(title="Detected genes per sample", x="", y="Genes with count > 0")
        + theme.GGPLOT_THEME
    )
    _rotate_x(p, counts.shape[1]).save(outdir / "detected_genes.png", verbose=False)

    # --- count-distribution densities (log scale) --------------------------
    # Subsample genes so the melt stays cheap on large cohorts (781 samples).
    max_genes = 3000
    sub = counts
    if counts.shape[0] > max_genes:
        sub = counts.sample(max_genes, random_state=config.RANDOM_STATE)
    long = (
        sub.replace(0, np.nan)
        .reset_index()
        .melt(id_vars="gene", var_name="sample", value_name="count")
        .dropna()
    )
    long = long.merge(
        meta["condition"].astype(str).rename("condition"),
        left_on="sample", right_index=True,
    )
    # colour by condition (not by 781 individual samples) and group by sample
    p = (
        ggplot(long, aes("count", color="condition", group="sample"))
        + geom_density(size=0.3, alpha=0.5)
        + scale_x_log10()
        + scale_color_manual(values=pal)
        + labs(title="Per-sample count distributions", x="Count (log10)", y="Density")
        + theme.GGPLOT_THEME
    )
    p.save(outdir / "count_densities.png", verbose=False)
    print(f"[p1] QC figures -> {outdir}")


def _rotate_x(p, n_samples: int | None = None):
    from plotnine import element_blank, element_text, theme as pn_theme
    # With many samples the x labels are unreadable; hide them.
    if n_samples is not None and n_samples > 40:
        return p + pn_theme(axis_text_x=element_blank())
    return p + pn_theme(axis_text_x=element_text(rotation=90, ha="center", size=7))


# --------------------------------------------------------------------------- #
def main() -> None:
    spec = config.design()
    counts = io_utils.load_counts("counts_raw.parquet")
    meta = io_utils.load_metadata("metadata_raw.parquet")

    meta = _tidy_condition(meta, spec)
    meta = _tidy_time(meta, spec)
    counts, meta = io_utils.align_counts_metadata(counts, meta)
    counts = filter_low_counts(counts)

    qc_plots(counts, meta)

    io_utils.save_counts(counts, "counts.parquet")
    io_utils.save_metadata(meta, "metadata.parquet")
    print("[p1] wrangling complete.")
    print(meta[["condition"]].value_counts())


if __name__ == "__main__":
    main()
