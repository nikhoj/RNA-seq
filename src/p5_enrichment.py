"""
Phase 5 - Functional / pathway analysis.

- Over-representation analysis (ORA) on DE gene sets via Enrichr (gseapy).
- GSEA (pre-ranked) on the full log2FC-ranked gene list.
- If temporal clusters exist (Phase 4), ORA is run per cluster to give each
  trajectory a biological identity.

Visuals: enrichment dot/bar plots, GSEA running-enrichment curves.

Requires network access for Enrichr / GSEA gene-set libraries.

Run:  RNASEQ_DATASET=airway python -m src.p5_enrichment
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402
from src import io_utils, theme  # noqa: E402,F401

GENE_SETS = ["GO_Biological_Process_2021", "KEGG_2021_Human", "MSigDB_Hallmark_2020"]


# --------------------------------------------------------------------------- #
def _load_de_tables() -> dict[str, pd.DataFrame]:
    """Load every DE results CSV produced by Phase 3 for this dataset."""
    tables = {}
    for p in sorted(config.TABLES.glob(f"{config.DATASET}_DE_*.csv")):
        name = p.stem.replace(f"{config.DATASET}_DE_", "")
        tables[name] = pd.read_csv(p, index_col=0)
    return tables


def ora(symbols: list[str], tag: str) -> pd.DataFrame | None:
    """Run Enrichr ORA for a gene-symbol list; save + plot the top terms."""
    import gseapy as gp

    symbols = [s for s in symbols if isinstance(s, str) and s and not s.isdigit()]
    if len(symbols) < 5:
        print(f"[p5] ORA {tag}: too few mapped symbols ({len(symbols)}) -- skipped")
        return None
    try:
        enr = gp.enrichr(gene_list=symbols, gene_sets=GENE_SETS, outdir=None)
    except Exception as exc:  # noqa: BLE001
        print(f"[p5] ORA {tag} failed: {exc}")
        return None
    res = enr.results.sort_values("Adjusted P-value")
    io_utils.save_table(res, f"ORA_{tag}.csv", index=False)
    _dotplot(res, tag)
    print(f"[p5] ORA {tag}: {len(res)} terms; top = {res.iloc[0]['Term'][:50]}")
    return res


def _dotplot(res: pd.DataFrame, tag: str, top: int = 12) -> None:
    from plotnine import (
        aes, geom_point, ggplot, labs, scale_color_cmap, scale_size,
        element_text, theme as pn_theme,
    )

    outdir = config.fig_dir("phase5_enrichment")
    d = res.head(top).copy()
    d["neglog10padj"] = -np.log10(d["Adjusted P-value"].clip(lower=1e-300))
    d["gene_count"] = d["Overlap"].str.split("/").str[0].astype(int)
    d["Term"] = d["Term"].str.slice(0, 55)
    d["Term"] = pd.Categorical(d["Term"], categories=d["Term"][::-1])

    p = (
        ggplot(d, aes("neglog10padj", "Term", size="gene_count", color="neglog10padj"))
        + geom_point()
        + scale_color_cmap(cmap_name="viridis")
        + scale_size(range=(2, 8))
        + labs(title=f"Enrichment (ORA): {tag}", x="-log10 adjusted p", y="")
        + theme.GGPLOT_THEME
        + pn_theme(axis_text_y=element_text(size=7))
    )
    p.save(outdir / f"ora_dotplot_{tag}.png", verbose=False)
    print(f"[p5] ORA dotplot -> {outdir/('ora_dotplot_'+tag+'.png')}")


def gsea_prerank(res: pd.DataFrame, tag: str) -> None:
    """Pre-ranked GSEA on the log2FC-ranked list; plot top pathway curves."""
    import gseapy as gp

    d = res.dropna(subset=["log2FoldChange", "Symbol"]).copy()
    d = d[~d["Symbol"].astype(str).str.isdigit()]
    rank = (
        d.groupby("Symbol")["log2FoldChange"].mean()
        .sort_values(ascending=False)
    )
    try:
        pre = gp.prerank(
            rnk=rank, gene_sets="MSigDB_Hallmark_2020",
            min_size=5, max_size=500, permutation_num=100,
            seed=config.RANDOM_STATE, outdir=None, verbose=False,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"[p5] GSEA {tag} failed: {exc}")
        return
    out = pre.res2d.sort_values("NES", key=lambda s: s.abs(), ascending=False)
    io_utils.save_table(out, f"GSEA_{tag}.csv", index=False)
    print(f"[p5] GSEA {tag}: top pathway = {out.iloc[0]['Term']} (NES={out.iloc[0]['NES']})")

    # running-enrichment curves for the two strongest terms
    outdir = config.fig_dir("phase5_enrichment")
    import matplotlib.pyplot as plt
    for term in out["Term"].head(2):
        try:
            ax = pre.plot(terms=term)
            fig = ax if hasattr(ax, "savefig") else plt.gcf()
            safe = term.replace(" ", "_").replace("/", "-")[:40]
            fig.savefig(outdir / f"gsea_{tag}_{safe}.png", dpi=150, bbox_inches="tight")
            plt.close("all")
        except Exception as exc:  # noqa: BLE001
            print(f"[p5] GSEA plot for {term} skipped: {exc}")
    print(f"[p5] GSEA curves -> {outdir}")


def per_cluster_ora() -> None:
    """Run ORA for each temporal cluster (Phase 4 output), if present."""
    cpath = config.TABLES / f"{config.DATASET}_temporal_clusters.csv"
    apath = config.processed_path("annotation.parquet")
    if not cpath.exists():
        print("[p5] no temporal clusters found -- skipping per-cluster ORA "
              "(run Phase 4 first for a time-course dataset).")
        return
    clusters = pd.read_csv(cpath, index_col=0)
    annot = pd.read_parquet(apath) if apath.exists() else None
    for cl in sorted(clusters["cluster"].unique()):
        gene_ids = clusters.index[clusters["cluster"] == cl].astype(str)
        if annot is not None:
            sym = annot["Symbol"]
            sym.index = sym.index.astype(str)
            symbols = gene_ids.map(sym).dropna().tolist()
        else:
            symbols = gene_ids.tolist()
        ora(symbols, f"cluster{cl}")


# --------------------------------------------------------------------------- #
def main():
    tables = _load_de_tables()
    if not tables:
        print("[p5] no DE tables found -- run Phase 3 first.")
        return

    for name, res in tables.items():
        sig = res[
            (res["padj"] < config.PADJ_THRESHOLD)
            & (res["log2FoldChange"].abs() >= config.LFC_THRESHOLD)
        ]
        symbols = sig["Symbol"].dropna().astype(str).tolist()
        ora(symbols, name)
        gsea_prerank(res, name)

    per_cluster_ora()
    print("[p5] functional enrichment complete.")


if __name__ == "__main__":
    main()
