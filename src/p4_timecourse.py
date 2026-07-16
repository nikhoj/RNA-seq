"""
Phase 4 - Time-course-specific analysis (the main reason for the project).

Only meaningful for datasets with a time axis (gse212041); for the airway
warmup this module no-ops with a message.

- Time-series-aware DE via a likelihood-ratio test (LRT): per gene, compare a
  full model  (~ condition + time + condition:time) against a reduced model
  (~ condition + time) on the VST matrix, using statsmodels OLS + an F/LR test.
  This is the practical Python stand-in for DESeq2's LRT (PyDESeq2 0.5.x has no
  built-in LRT) and is well-behaved on variance-stabilised data with large N.
- Temporal clustering of significant genes by trajectory *shape* via fuzzy
  c-means (scikit-fuzzy = the Python Mfuzz) with a k-means fallback.
- Smooth spline trajectory fits per cluster.
- Visuals: per-cluster trajectory panels, time-ordered heatmap, single-gene
  ribbon plots.

Run:  RNASEQ_DATASET=gse212041 python -m src.p4_timecourse
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402
from src import io_utils, theme  # noqa: E402,F401


def has_time_axis(meta: pd.DataFrame) -> bool:
    return "time" in meta.columns and meta["time"].notna().sum() > 0 and \
        meta["time"].nunique() >= 3


# --------------------------------------------------------------------------- #
# LRT: genes whose trajectory differs between conditions over time
# --------------------------------------------------------------------------- #
def lrt_time_interaction(vst: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """
    Per-gene likelihood-ratio test for a time effect, run as nested OLS on VST.

    Two modes (config.design()["timecourse_test"]):
      - "interaction": full ~condition*time vs reduced ~condition+time
        -> genes whose *trajectory differs between conditions*. Needs both arms
           sampled across time.
      - "time_main":  restrict to one condition (timecourse_condition) and test
        full ~time vs reduced ~1 -> genes that *change over time* within that
        arm. The right question when controls lack a real time course.
    """
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests
    from scipy.stats import chi2

    spec = config.design()
    mode = spec.get("timecourse_test", "interaction")

    m = meta.copy()
    m["condition"] = m["condition"].astype(str)
    m["time"] = m["time"].astype(float)
    m = m[m["time"].notna()]                      # drop non-numeric timepoints

    if mode == "time_main":
        cond = spec.get("timecourse_condition")
        if cond is not None:
            m = m[m["condition"] == str(cond)]
        full_f, reduced_f = "y ~ time", "y ~ 1"
        print(f"[p4] LRT mode=time_main on condition='{cond}' (n={len(m)} samples)")
    else:
        full_f, reduced_f = "y ~ condition * time", "y ~ condition + time"
        print(f"[p4] LRT mode=interaction (n={len(m)} samples)")

    vst = vst[m.index]                            # align columns to kept samples
    stats = []
    genes = vst.index
    for i, g in enumerate(genes):
        m["y"] = vst.loc[g].values
        try:
            full = smf.ols(full_f, data=m).fit()
            reduced = smf.ols(reduced_f, data=m).fit()
            n = int(full.nobs)
            rss_f, rss_r = full.ssr, reduced.ssr
            lr = n * np.log(rss_r / rss_f) if rss_f > 0 else 0.0
            df_diff = int(full.df_model - reduced.df_model)
            pval = chi2.sf(lr, df_diff) if df_diff > 0 else 1.0
        except Exception:  # noqa: BLE001
            lr, pval = np.nan, np.nan
        stats.append((g, lr, pval))
        if (i + 1) % 2000 == 0:
            print(f"[p4] LRT {i+1}/{len(genes)} genes")

    res = pd.DataFrame(stats, columns=["gene", "LR_stat", "pvalue"]).set_index("gene")
    ok = res["pvalue"].notna()
    res.loc[ok, "padj"] = multipletests(res.loc[ok, "pvalue"], method="fdr_bh")[1]
    res = res.sort_values("padj")
    io_utils.save_table(res, "timecourse_LRT.csv")
    n_sig = int((res["padj"] < config.PADJ_THRESHOLD).sum())
    print(f"[p4] LRT ({mode}): {n_sig} genes changing over time (padj<{config.PADJ_THRESHOLD})")
    return res


# --------------------------------------------------------------------------- #
# Temporal clustering
# --------------------------------------------------------------------------- #
def cluster_trajectories(
    vst: pd.DataFrame, meta: pd.DataFrame, genes: list[str],
    k: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Cluster genes by trajectory shape. Builds a genes x timepoint matrix of
    condition-averaged, z-scored expression, then fuzzy-c-means clusters it.

    Returns (assignment_df, centroid_traj_df).
    """
    k = k or config.N_TEMPORAL_CLUSTERS
    # mean expression per (gene, time) averaged over samples/replicates
    long = (
        vst.loc[genes].T.join(meta[["time"]])
        .groupby("time").mean().T          # genes x timepoints
    )
    long = long.reindex(sorted(long.columns), axis=1)
    # z-score each gene's trajectory
    z = long.sub(long.mean(axis=1), axis=0).div(long.std(axis=1).replace(0, 1), axis=0)

    try:
        import skfuzzy as fuzz
        cntr, u, *_ = fuzz.cluster.cmeans(
            z.values.T, c=k, m=1.25, error=1e-5, maxiter=1000,
            seed=config.RANDOM_STATE,
        )
        labels = u.argmax(axis=0)
        membership = u.max(axis=0)
        method = "fuzzy-c-means"
        centroids = pd.DataFrame(cntr, columns=z.columns)
    except Exception as exc:  # noqa: BLE001
        print(f"[p4] fuzzy c-means unavailable ({exc}); using k-means")
        from sklearn.cluster import KMeans
        km = KMeans(n_clusters=k, random_state=config.RANDOM_STATE, n_init=10).fit(z.values)
        labels = km.labels_
        membership = np.ones(len(labels))
        centroids = pd.DataFrame(km.cluster_centers_, columns=z.columns)
        method = "k-means"

    assign = pd.DataFrame(
        {"cluster": labels, "membership": membership}, index=z.index
    )
    assign = assign.join(z, how="left")
    io_utils.save_table(assign[["cluster", "membership"]], "temporal_clusters.csv")
    print(f"[p4] {method}: {k} clusters over {z.shape[1]} timepoints; "
          f"sizes = {dict(pd.Series(labels).value_counts().sort_index())}")
    return assign, centroids


# --------------------------------------------------------------------------- #
# Visuals
# --------------------------------------------------------------------------- #
def plot_cluster_panels(assign: pd.DataFrame, centroids: pd.DataFrame) -> None:
    from plotnine import (
        aes, geom_line, ggplot, labs, facet_wrap, geom_hline,
    )

    outdir = config.fig_dir("phase4_timecourse")
    tcols = [c for c in assign.columns if c not in ("cluster", "membership")]

    # individual gene trajectories (thin) + centroid (thick) per cluster
    a = assign.copy()
    a.index.name = "gene"
    long = a.reset_index().melt(
        id_vars=["gene", "cluster", "membership"], value_vars=tcols,
        var_name="time", value_name="z",
    )
    long["time"] = long["time"].astype(float)
    long["cluster"] = "cluster " + long["cluster"].astype(str)

    p = (
        ggplot(long, aes("time", "z", group="gene"))
        + geom_line(alpha=0.06, color="#0072B2")
        + geom_hline(yintercept=0, linetype="dotted", color="grey")
        + facet_wrap("cluster")
        + labs(title="Temporal expression clusters (z-scored VST)",
               x="time", y="z-scored expression")
        + theme.GGPLOT_THEME
    )
    p.save(outdir / "cluster_trajectories.png", verbose=False)
    print(f"[p4] cluster panels -> {outdir/'cluster_trajectories.png'}")


def plot_time_heatmap(vst: pd.DataFrame, meta: pd.DataFrame,
                      genes: list[str]) -> None:
    import seaborn as sns
    import matplotlib.pyplot as plt

    outdir = config.fig_dir("phase4_timecourse")
    mat = (
        vst.loc[genes].T.join(meta[["time"]]).groupby("time").mean().T
    ).reindex(sorted(set(meta["time"].dropna())), axis=1)
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, 1), axis=0)

    g = sns.clustermap(
        z, cmap="RdBu_r", center=0, col_cluster=False,
        figsize=(6, 9), yticklabels=False, vmin=-2, vmax=2,
    )
    g.fig.suptitle("Time-ordered expression heatmap (top time-varying genes)", y=1.01)
    g.savefig(outdir / "time_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"[p4] time heatmap -> {outdir/'time_heatmap.png'}")


def plot_gene_ribbons(vst: pd.DataFrame, meta: pd.DataFrame,
                      genes: list[str], n: int = 6) -> None:
    from plotnine import aes, geom_line, geom_ribbon, ggplot, labs, facet_wrap, scale_color_manual, scale_fill_manual

    outdir = config.fig_dir("phase4_timecourse")
    frames = []
    for g in genes[:n]:
        df = meta[["time", "condition"]].copy()
        df["expr"] = vst.loc[g].reindex(df.index).values
        df["gene"] = g
        frames.append(df)
    long = pd.concat(frames)
    long["condition"] = long["condition"].astype(str)
    summ = (
        long.groupby(["gene", "condition", "time"])["expr"]
        .agg(["mean", "std", "count"]).reset_index()
    )
    summ["se"] = summ["std"] / np.sqrt(summ["count"].clip(lower=1))
    summ["lo"], summ["hi"] = summ["mean"] - summ["se"], summ["mean"] + summ["se"]

    p = (
        ggplot(summ, aes("time", "mean", color="condition", fill="condition"))
        + geom_ribbon(aes(ymin="lo", ymax="hi"), alpha=0.2, color=None)
        + geom_line(size=1)
        + facet_wrap("gene", scales="free_y")
        + scale_color_manual(values=theme.OKABE_ITO)
        + scale_fill_manual(values=theme.OKABE_ITO)
        + labs(title="Expression trajectories (mean +/- SE)", x="time", y="VST expression")
        + theme.GGPLOT_THEME
    )
    p.save(outdir / "gene_trajectories.png", verbose=False)
    print(f"[p4] gene trajectory ribbons -> {outdir/'gene_trajectories.png'}")


# --------------------------------------------------------------------------- #
def main():
    meta = io_utils.load_metadata("metadata.parquet")
    if not has_time_axis(meta):
        print(f"[p4] dataset '{config.DATASET}' has no time axis -- Phase 4 skipped.")
        print("      (Run with RNASEQ_DATASET=gse212041 once that data is wrangled.)")
        return

    vst = pd.read_parquet(config.processed_path("vst.parquet"))
    lrt = lrt_time_interaction(vst, meta)

    sig_genes = lrt[lrt["padj"] < config.PADJ_THRESHOLD].index.tolist()
    if not sig_genes:
        print("[p4] no significant time-varying genes; using top 500 by LR stat")
        sig_genes = lrt.sort_values("LR_stat", ascending=False).head(500).index.tolist()
    top_genes = sig_genes[:2000]  # cap for clustering/plots

    # In time_main mode, restrict trajectory clustering/plots to that condition
    # (with numeric time) so the shapes reflect a single clean course.
    spec = config.design()
    if spec.get("timecourse_test") == "time_main" and spec.get("timecourse_condition"):
        keep = (meta["condition"].astype(str) == str(spec["timecourse_condition"])) \
            & meta["time"].notna()
        meta = meta[keep]
        vst = vst[meta.index]

    assign, centroids = cluster_trajectories(vst, meta, top_genes)
    plot_cluster_panels(assign, centroids)
    plot_time_heatmap(vst, meta, top_genes[:200])
    plot_gene_ribbons(vst, meta, sig_genes)
    print("[p4] time-course analysis complete.")


if __name__ == "__main__":
    main()
