"""
Phase 3 - Differential expression (core).

- fit the DESeq2 model and run Wald contrasts (each non-baseline condition level
  vs the baseline)
- LFC shrinkage (apeGLM-style via PyDESeq2)
- significance calls (padj < PADJ_THRESHOLD & |log2FC| >= LFC_THRESHOLD)
- visuals: volcano (plotnine + adjustText), MA plot, top-gene boxplots,
  up/down bar summary

Run:  RNASEQ_DATASET=airway python -m src.p3_de
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402
from src import io_utils, theme  # noqa: E402,F401


def _annotate_symbols(df: pd.DataFrame) -> pd.DataFrame:
    """Attach gene Symbol from the cached annotation, if available."""
    apath = config.processed_path("annotation.parquet")
    if apath.exists():
        annot = pd.read_parquet(apath)
        # anndata stringifies gene IDs; annotation index is int. Map on string.
        sym = annot["Symbol"]
        sym.index = sym.index.astype(str)
        df["Symbol"] = df.index.astype(str).map(sym)
        df["Symbol"] = df["Symbol"].fillna(df.index.to_series().astype(str))
    else:
        df["Symbol"] = df.index.astype(str)
    return df


def run_contrasts(dds) -> dict[str, pd.DataFrame]:
    """
    Run one Wald contrast per non-baseline condition level vs baseline.
    Returns {contrast_name: results_dataframe}.
    """
    from pydeseq2.ds import DeseqStats

    levels = list(dds.obs["condition"].cat.categories)
    baseline = levels[0]
    results: dict[str, pd.DataFrame] = {}

    for level in levels[1:]:
        contrast = ["condition", level, baseline]
        name = f"{level}_vs_{baseline}"
        stat = DeseqStats(dds, contrast=contrast, quiet=True)
        stat.summary()
        try:
            stat.lfc_shrink(coeff=f"condition[T.{level}]")
        except Exception as exc:  # noqa: BLE001
            print(f"[p3] lfc_shrink skipped for {name}: {exc}")
        res = stat.results_df.copy()
        res = _annotate_symbols(res)
        res["significant"] = (
            (res["padj"] < config.PADJ_THRESHOLD)
            & (res["log2FoldChange"].abs() >= config.LFC_THRESHOLD)
        )
        res["direction"] = np.where(res["log2FoldChange"] >= 0, "up", "down")
        results[name] = res
        io_utils.save_table(res, f"DE_{name}.csv")
        n_sig = int(res["significant"].sum())
        print(f"[p3] {name}: {n_sig} significant genes "
              f"(padj<{config.PADJ_THRESHOLD}, |LFC|>={config.LFC_THRESHOLD})")
    return results


# --------------------------------------------------------------------------- #
# Visuals
# --------------------------------------------------------------------------- #
def volcano(res: pd.DataFrame, name: str, n_labels: int = 12) -> None:
    from plotnine import (
        aes, geom_point, geom_hline, geom_vline, ggplot, labs,
        scale_color_manual,
    )
    import matplotlib.pyplot as plt
    from adjustText import adjust_text

    outdir = config.fig_dir("phase3_de")
    d = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    d["neglog10padj"] = -np.log10(d["padj"].clip(lower=1e-300))
    d["status"] = np.where(
        ~d["significant"], "ns",
        np.where(d["log2FoldChange"] > 0, "up", "down"),
    )
    palette = {"ns": "#b0b0b0", "up": "#D55E00", "down": "#0072B2"}

    p = (
        ggplot(d, aes("log2FoldChange", "neglog10padj", color="status"))
        + geom_point(size=1.2, alpha=0.6)
        + geom_hline(yintercept=-np.log10(config.PADJ_THRESHOLD),
                     linetype="dashed", color="grey")
        + geom_vline(xintercept=[-config.LFC_THRESHOLD, config.LFC_THRESHOLD],
                     linetype="dashed", color="grey")
        + scale_color_manual(values=palette)
        + labs(title=f"Volcano: {name}", x="log2 fold change",
               y="-log10 adjusted p")
        + theme.GGPLOT_THEME
    )
    fig = p.draw()
    ax = fig.axes[0]
    top = d[d["significant"]].reindex(
        d[d["significant"]]["neglog10padj"].sort_values(ascending=False).index
    ).head(n_labels)
    texts = [
        ax.text(r["log2FoldChange"], r["neglog10padj"], r["Symbol"], fontsize=7)
        for _, r in top.iterrows()
    ]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.4))
    fig.savefig(outdir / f"volcano_{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[p3] volcano -> {outdir/('volcano_'+name+'.png')}")


def ma_plot(res: pd.DataFrame, name: str) -> None:
    from plotnine import aes, geom_point, geom_hline, ggplot, labs, scale_x_log10, scale_color_manual

    outdir = config.fig_dir("phase3_de")
    d = res.dropna(subset=["baseMean", "log2FoldChange"]).copy()
    d = d[d["baseMean"] > 0]
    d["sig"] = np.where(d["significant"], "significant", "ns")
    p = (
        ggplot(d, aes("baseMean", "log2FoldChange", color="sig"))
        + geom_point(size=0.8, alpha=0.5)
        + geom_hline(yintercept=0, color="black", size=0.4)
        + scale_x_log10()
        + scale_color_manual(values={"significant": "#D55E00", "ns": "#b0b0b0"})
        + labs(title=f"MA plot: {name}", x="mean of normalized counts",
               y="log2 fold change")
        + theme.GGPLOT_THEME
    )
    p.save(outdir / f"ma_{name}.png", verbose=False)
    print(f"[p3] MA -> {outdir/('ma_'+name+'.png')}")


def top_gene_boxplots(dds, res: pd.DataFrame, meta: pd.DataFrame, name: str,
                      n: int = 6) -> None:
    from plotnine import (
        aes, geom_boxplot, geom_jitter, ggplot, labs, facet_wrap,
        scale_fill_manual,
    )

    outdir = config.fig_dir("phase3_de")
    norm = pd.DataFrame(
        dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names
    )  # samples x genes
    top = res[res["significant"]].reindex(
        res[res["significant"]]["padj"].sort_values().index
    ).head(n)
    if top.empty:
        print(f"[p3] no significant genes to box-plot for {name}")
        return

    frames = []
    for gid, r in top.iterrows():
        vals = np.log2(norm[gid] + 1)
        frames.append(pd.DataFrame({
            "log2_norm": vals.values,
            "condition": meta["condition"].astype(str).values,
            "gene": f"{r['Symbol']}",
        }))
    long = pd.concat(frames)
    long["gene"] = pd.Categorical(long["gene"], categories=top["Symbol"].tolist())

    p = (
        ggplot(long, aes("condition", "log2_norm", fill="condition"))
        + geom_boxplot(outlier_alpha=0)
        + geom_jitter(width=0.15, size=1.5, alpha=0.7)
        + facet_wrap("gene", scales="free_y")
        + scale_fill_manual(values=theme.OKABE_ITO)
        + labs(title=f"Top DE genes: {name}", x="", y="log2(normalized count + 1)")
        + theme.GGPLOT_THEME
    )
    p.save(outdir / f"top_genes_{name}.png", verbose=False)
    print(f"[p3] top-gene boxplots -> {outdir/('top_genes_'+name+'.png')}")


def de_summary_bar(results: dict[str, pd.DataFrame]) -> None:
    from plotnine import aes, geom_col, geom_text, ggplot, labs, position_dodge, scale_fill_manual

    outdir = config.fig_dir("phase3_de")
    rows = []
    for name, res in results.items():
        sig = res[res["significant"]]
        rows.append({"contrast": name, "direction": "up",
                     "n": int((sig["direction"] == "up").sum())})
        rows.append({"contrast": name, "direction": "down",
                     "n": int((sig["direction"] == "down").sum())})
    df = pd.DataFrame(rows)
    p = (
        ggplot(df, aes("contrast", "n", fill="direction"))
        + geom_col(position=position_dodge(width=0.8), width=0.7)
        + geom_text(aes(label="n"), position=position_dodge(width=0.8),
                    va="bottom", size=8)
        + scale_fill_manual(values={"up": "#D55E00", "down": "#0072B2"})
        + labs(title="DE gene counts per contrast", x="", y="# significant genes")
        + theme.GGPLOT_THEME
    )
    p.save(outdir / "de_summary_bar.png", verbose=False)
    print(f"[p3] DE summary bar -> {outdir/'de_summary_bar.png'}")


# --------------------------------------------------------------------------- #
def main():
    from src.p2_qc import build_dds

    counts = io_utils.load_counts("counts.parquet")
    meta = io_utils.load_metadata("metadata.parquet")
    dds = build_dds(counts, meta)

    results = run_contrasts(dds)
    for name, res in results.items():
        volcano(res, name)
        ma_plot(res, name)
        top_gene_boxplots(dds, res, meta, name)
    de_summary_bar(results)
    print("[p3] differential expression complete.")
    return dds, results


if __name__ == "__main__":
    main()
