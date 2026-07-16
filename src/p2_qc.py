"""
Phase 2 - Exploratory data analysis (QC).

- DESeq2 median-of-ratios size factors (via PyDESeq2)
- variance-stabilising transform (VST) for visualisation
- PCA on the VST matrix
- sample-to-sample distance clustermap
- normalized-count boxplots (before vs after)

Run:  RNASEQ_DATASET=airway python -m src.p2_qc
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402
from src import io_utils, theme  # noqa: E402,F401


def build_dds(counts: pd.DataFrame, meta: pd.DataFrame):
    """
    Build and fit a PyDESeq2 DeseqDataSet (size factors + dispersions + VST).
    PyDESeq2 wants samples x genes counts and integer values.
    """
    from pydeseq2.dds import DeseqDataSet

    counts_sxg = counts.T.astype(int)  # samples x genes
    factors = [f for f in config.design()["design_factors"] if f in meta.columns]
    design = "~ " + " + ".join(factors)

    dds = DeseqDataSet(
        counts=counts_sxg,
        metadata=meta,
        design=design,
        refit_cooks=True,
        quiet=True,
    )
    dds.deseq2()          # size factors, dispersions, LFCs
    dds.vst(use_design=False)
    print(f"[p2] DeseqDataSet fitted; design = {design}")
    return dds


def vst_matrix(dds) -> pd.DataFrame:
    """Return the VST matrix as genes x samples."""
    vst = pd.DataFrame(
        dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names
    ).T  # genes x samples
    return vst


# --------------------------------------------------------------------------- #
# Visuals
# --------------------------------------------------------------------------- #
def plot_pca(vst: pd.DataFrame, meta: pd.DataFrame) -> None:
    from sklearn.decomposition import PCA
    from plotnine import aes, geom_point, ggplot, labs, scale_color_manual

    outdir = config.fig_dir("phase2_qc")
    top = vst.var(axis=1).sort_values(ascending=False).head(500).index
    X = vst.loc[top].T.values  # samples x top-genes
    pca = PCA(n_components=2, random_state=config.RANDOM_STATE).fit(X)
    pcs = pca.transform(X)
    var = pca.explained_variance_ratio_ * 100

    df = meta.copy()
    df["PC1"], df["PC2"] = pcs[:, 0], pcs[:, 1]
    df["condition"] = df["condition"].astype(str)
    has_time = "time" in df.columns and df["time"].notna().any()

    aesthetic = aes("PC1", "PC2", color="condition")
    if has_time:
        df["time_factor"] = df["time_factor"].astype(str)
        aesthetic = aes("PC1", "PC2", color="time_factor", shape="condition")

    p = (
        ggplot(df, aesthetic)
        + geom_point(size=4, alpha=0.85)
        + scale_color_manual(values=theme.OKABE_ITO)
        + labs(
            title="PCA on VST (top 500 variable genes)",
            x=f"PC1 ({var[0]:.1f}%)", y=f"PC2 ({var[1]:.1f}%)",
        )
        + theme.GGPLOT_THEME
    )
    p.save(outdir / "pca.png", verbose=False)
    print(f"[p2] PCA -> {outdir/'pca.png'}")


def plot_sample_distance(vst: pd.DataFrame, meta: pd.DataFrame) -> None:
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.spatial.distance import pdist, squareform

    outdir = config.fig_dir("phase2_qc")
    d = squareform(pdist(vst.T.values, metric="euclidean"))
    dist = pd.DataFrame(d, index=vst.columns, columns=vst.columns)

    labels = meta["condition"].astype(str)
    lut = dict(zip(labels.unique(), theme.OKABE_ITO))
    colors = labels.map(lut)

    show_labels = dist.shape[0] <= 40  # unreadable for large cohorts
    g = sns.clustermap(
        dist, cmap="mako_r", figsize=(8, 8),
        row_colors=colors.values, col_colors=colors.values,
        xticklabels=show_labels, yticklabels=show_labels,
    )
    g.fig.suptitle("Sample-to-sample distances (VST, Euclidean)", y=1.02)
    g.savefig(outdir / "sample_distance_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"[p2] distance heatmap -> {outdir/'sample_distance_heatmap.png'}")


def plot_normalization_boxplots(counts: pd.DataFrame, dds, meta: pd.DataFrame) -> None:
    from plotnine import (
        aes, geom_boxplot, ggplot, labs, facet_wrap, element_blank,
        element_text, theme as pn_theme,
    )

    outdir = config.fig_dir("phase2_qc")
    norm = pd.DataFrame(
        dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names
    ).T  # genes x samples

    # Subsample genes so the melt stays cheap on large cohorts.
    max_genes = 3000
    if counts.shape[0] > max_genes:
        genes = counts.sample(max_genes, random_state=config.RANDOM_STATE).index
        counts, norm = counts.loc[genes], norm.loc[genes]

    def _long(mat, tag):
        m = np.log2(mat.replace(0, np.nan) + 1)
        long = m.reset_index().melt(id_vars="gene", var_name="sample",
                                    value_name="log2_count").dropna()
        long["stage"] = tag
        return long

    long = pd.concat([_long(counts, "raw"), _long(norm, "normalized")])
    long["condition"] = long["sample"].map(meta["condition"].astype(str))

    p = (
        ggplot(long, aes("sample", "log2_count", fill="condition"))
        + geom_boxplot(outlier_size=0.2, outlier_alpha=0.15)
        + facet_wrap("stage", ncol=1)
        + labs(title="Per-sample counts: raw vs normalized",
               x="", y="log2(count + 1)")
        + theme.GGPLOT_THEME
        + pn_theme(axis_text_x=(
            element_blank() if norm.shape[1] > 40
            else element_text(rotation=90, size=7)))
    )
    p.save(outdir / "normalization_boxplots.png", verbose=False)
    print(f"[p2] normalization boxplots -> {outdir/'normalization_boxplots.png'}")


# --------------------------------------------------------------------------- #
def main():
    counts = io_utils.load_counts("counts.parquet")
    meta = io_utils.load_metadata("metadata.parquet")

    dds = build_dds(counts, meta)
    vst = vst_matrix(dds)

    # cache VST + size factors for later phases
    vst.to_parquet(config.processed_path("vst.parquet"))
    dds.obs["size_factors"].rename("size_factor").to_frame().to_parquet(
        config.processed_path("size_factors.parquet")
    )

    plot_pca(vst, meta)
    plot_sample_distance(vst, meta)
    plot_normalization_boxplots(counts, dds, meta)
    print("[p2] EDA/QC complete.")
    return dds


if __name__ == "__main__":
    main()
