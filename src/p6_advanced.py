"""
Phase 6 - Advanced & integrative (optional).

- Gene co-expression modules (WGCNA-style): correlation network on the most
  variable genes, hierarchical module detection, module-eigengene summary.
- Signature scoring over time: score a gene set per sample (mean z-score) and
  plot the score across time / condition (needs a time axis).
- Optional ML: classify a binary outcome from expression with logistic
  regression + ROC (needs an outcome label in the metadata).

Each block guards for its prerequisites, so it degrades gracefully on the
airway warmup (no time, no outcome).

Run:  RNASEQ_DATASET=gse212041 python -m src.p6_advanced
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
# Co-expression modules
# --------------------------------------------------------------------------- #
def coexpression_modules(vst: pd.DataFrame, n_genes: int = 1500,
                         n_modules: int = 8) -> pd.DataFrame:
    """
    WGCNA-style modules: correlation-distance hierarchical clustering of the top
    variable genes into modules. Returns a gene->module assignment.
    """
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    top = vst.var(axis=1).sort_values(ascending=False).head(n_genes).index
    X = vst.loc[top]
    corr = np.corrcoef(X.values)                 # gene x gene
    dist = 1 - corr
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2                    # enforce symmetry
    Z = linkage(squareform(dist, checks=False), method="average")
    modules = fcluster(Z, t=n_modules, criterion="maxclust")

    assign = pd.DataFrame({"module": modules}, index=top)
    io_utils.save_table(assign, "coexpression_modules.csv")
    print(f"[p6] co-expression: {n_modules} modules over {n_genes} genes; "
          f"sizes = {dict(pd.Series(modules).value_counts().sort_index())}")
    _plot_module_heatmap(X, assign)
    return assign


def _plot_module_heatmap(X: pd.DataFrame, assign: pd.DataFrame) -> None:
    import seaborn as sns
    import matplotlib.pyplot as plt

    outdir = config.fig_dir("phase4_timecourse")  # reuse advanced-figure folder
    order = assign.sort_values("module").index
    z = X.loc[order]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, 1), axis=0)
    lut = dict(zip(sorted(assign["module"].unique()), theme.OKABE_ITO * 3))
    row_colors = assign.loc[order, "module"].map(lut)
    g = sns.clustermap(
        z, cmap="RdBu_r", center=0, row_cluster=False, col_cluster=True,
        row_colors=row_colors.values, yticklabels=False, figsize=(7, 9),
        vmin=-2, vmax=2,
    )
    g.fig.suptitle("Co-expression modules (VST, z-scored)", y=1.01)
    g.savefig(outdir / "coexpression_modules.png", dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"[p6] module heatmap -> {outdir/'coexpression_modules.png'}")


# --------------------------------------------------------------------------- #
# Signature scoring over time
# --------------------------------------------------------------------------- #
def signature_over_time(vst: pd.DataFrame, meta: pd.DataFrame,
                        signature: list[str], name: str = "signature") -> None:
    if "time" not in meta.columns or meta["time"].notna().sum() == 0:
        print("[p6] signature-over-time skipped (no time axis).")
        return
    from plotnine import aes, geom_line, geom_ribbon, ggplot, labs, scale_color_manual, scale_fill_manual

    present = [g for g in signature if g in vst.index]
    if not present:
        print("[p6] none of the signature genes present -- skipped.")
        return
    sub = vst.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)
    score = z.mean(axis=0)                        # per-sample signature score

    df = meta[["time", "condition"]].copy()
    df["score"] = score.reindex(df.index).values
    df["condition"] = df["condition"].astype(str)
    summ = df.groupby(["condition", "time"])["score"].agg(["mean", "std", "count"]).reset_index()
    summ["se"] = summ["std"] / np.sqrt(summ["count"].clip(lower=1))
    summ["lo"], summ["hi"] = summ["mean"] - summ["se"], summ["mean"] + summ["se"]

    outdir = config.fig_dir("phase4_timecourse")
    p = (
        ggplot(summ, aes("time", "mean", color="condition", fill="condition"))
        + geom_ribbon(aes(ymin="lo", ymax="hi"), alpha=0.2, color=None)
        + geom_line(size=1.2)
        + scale_color_manual(values=theme.OKABE_ITO)
        + scale_fill_manual(values=theme.OKABE_ITO)
        + labs(title=f"{name} score over time", x="time", y="signature score (z)")
        + theme.GGPLOT_THEME
    )
    p.save(outdir / f"signature_{name}.png", verbose=False)
    print(f"[p6] signature-over-time -> {outdir/('signature_'+name+'.png')}")


# --------------------------------------------------------------------------- #
# Optional ML: outcome classification + ROC
# --------------------------------------------------------------------------- #
def classify_outcome(vst: pd.DataFrame, meta: pd.DataFrame,
                     outcome_col: str | None = None) -> None:
    spec = config.design()
    outcome_col = outcome_col or spec.get("outcome_col")
    if not outcome_col or outcome_col not in meta.columns:
        print(f"[p6] outcome column '{outcome_col}' absent -- ML skipped.")
        return
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_predict, StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.metrics import roc_curve, roc_auc_score
    from plotnine import aes, geom_line, geom_abline, ggplot, labs

    raw = meta[outcome_col].astype(str)
    # Binarise per config (e.g. acuity 1-2 = severe vs 4-5 = mild); otherwise
    # require the column to already be a 2-class label.
    pos, neg = spec.get("outcome_positive"), spec.get("outcome_negative")
    if pos and neg:
        mask = raw.isin([str(x) for x in pos] + [str(x) for x in neg])
        sub_meta = meta[mask]
        y_bin = raw[mask].isin([str(x) for x in pos]).astype(int).values
        label = f"{outcome_col}: {'/'.join(map(str,pos))} vs {'/'.join(map(str,neg))}"
    else:
        classes = raw.value_counts()
        if len(classes) != 2 or classes.min() < 5:
            print(f"[p6] outcome not a usable binary label ({dict(classes)}) -- skipped.")
            return
        sub_meta = meta
        y_bin = (raw == classes.index[0]).astype(int).values
        label = outcome_col

    if len(np.unique(y_bin)) < 2 or min(np.bincount(y_bin)) < 5:
        print(f"[p6] too few samples per outcome class -- ML skipped.")
        return
    top = vst.var(axis=1).sort_values(ascending=False).head(500).index
    X = vst.loc[top].T.reindex(sub_meta.index).values

    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, C=0.1))
    cv = StratifiedKFold(5, shuffle=True, random_state=config.RANDOM_STATE)
    proba = cross_val_predict(clf, X, y_bin, cv=cv, method="predict_proba")[:, 1]
    auc = roc_auc_score(y_bin, proba)
    fpr, tpr, _ = roc_curve(y_bin, proba)

    outdir = config.fig_dir("phase4_timecourse")
    roc = pd.DataFrame({"fpr": fpr, "tpr": tpr})
    p = (
        ggplot(roc, aes("fpr", "tpr"))
        + geom_abline(intercept=0, slope=1, linetype="dashed", color="grey")
        + geom_line(size=1.2, color="#D55E00")
        + labs(title=f"Outcome ROC ({label}), AUC={auc:.2f}",
               x="false positive rate", y="true positive rate")
        + theme.GGPLOT_THEME
    )
    p.save(outdir / "outcome_roc.png", verbose=False)
    print(f"[p6] outcome classifier ROC AUC={auc:.3f} -> {outdir/'outcome_roc.png'}")


# --------------------------------------------------------------------------- #
def main():
    vst = pd.read_parquet(config.processed_path("vst.parquet"))
    meta = io_utils.load_metadata("metadata.parquet")

    coexpression_modules(vst)

    # A small demo signature (glucocorticoid response) so the block runs on any
    # dataset that has the genes; swap for a domain-relevant set in real use.
    gc_signature = ["FKBP5", "TSC22D3", "ZBTB16", "PER1", "DUSP1", "KLF15"]
    # map symbols -> gene ids if the index is Entrez
    apath = config.processed_path("annotation.parquet")
    if apath.exists():
        annot = pd.read_parquet(apath)
        sym2id = {s: str(i) for i, s in annot["Symbol"].items()}
        gc_ids = [sym2id[s] for s in gc_signature if s in sym2id]
        gc_ids = [g for g in gc_ids if g in vst.index]
    else:
        gc_ids = [g for g in gc_signature if g in vst.index]
    signature_over_time(vst, meta, gc_ids, name="glucocorticoid_response")

    # outcome ML — uses config.design()["outcome_col"] (+ optional binarisation)
    classify_outcome(vst, meta)
    print("[p6] advanced analyses complete.")


if __name__ == "__main__":
    main()
