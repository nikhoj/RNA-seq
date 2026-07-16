"""
Phase 7 - Reporting.

Assembles the strongest figures and result tables into a single Markdown report
(report/<dataset>_report.md) with a per-phase narrative. Figures are referenced
by relative path so the report renders on GitHub or any Markdown viewer.

Run:  RNASEQ_DATASET=airway python -m src.p7_report
"""
from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402

PHASE_FIGURES = {
    "Phase 1 - QC": "phase1_qc",
    "Phase 2 - EDA / normalization": "phase2_qc",
    "Phase 3 - Differential expression": "phase3_de",
    "Phase 4 - Time-course": "phase4_timecourse",
    "Phase 5 - Enrichment": "phase5_enrichment",
}


def _rel(p: Path) -> str:
    return str(p.relative_to(config.ROOT))


def _de_highlights() -> str:
    lines = []
    for p in sorted(config.TABLES.glob(f"{config.DATASET}_DE_*.csv")):
        res = pd.read_csv(p, index_col=0)
        sig = res[(res["padj"] < config.PADJ_THRESHOLD)
                  & (res["log2FoldChange"].abs() >= config.LFC_THRESHOLD)]
        name = p.stem.replace(f"{config.DATASET}_DE_", "")
        up = int((sig["log2FoldChange"] > 0).sum())
        down = int((sig["log2FoldChange"] < 0).sum())
        top = res.dropna(subset=["padj"]).sort_values("padj").head(8)
        top_str = ", ".join(str(s) for s in top["Symbol"])
        lines.append(f"- **{name}**: {len(sig)} significant "
                     f"({up} up, {down} down). Top genes: {top_str}.")
    return "\n".join(lines) if lines else "_No DE tables found._"


def build_report() -> Path:
    config.REPORT.mkdir(parents=True, exist_ok=True)
    out = config.REPORT / f"{config.DATASET}_report.md"

    md: list[str] = []
    md.append(f"# RNA-seq analysis report - `{config.DATASET}`\n")
    md.append(f"_Generated {date.today().isoformat()} | "
              f"thresholds: padj < {config.PADJ_THRESHOLD}, "
              f"|log2FC| >= {config.LFC_THRESHOLD}_\n")

    md.append("## Differential expression highlights\n")
    md.append(_de_highlights() + "\n")

    for title, folder in PHASE_FIGURES.items():
        figdir = config.FIGURES / folder
        pngs = sorted(figdir.glob("*.png")) if figdir.exists() else []
        if not pngs:
            continue
        md.append(f"## {title}\n")
        for png in pngs:
            md.append(f"### {png.stem}\n")
            md.append(f"![{png.stem}]({_rel(png)})\n")

    md.append("## Tables\n")
    for csv in sorted(config.TABLES.glob(f"{config.DATASET}_*.csv")):
        md.append(f"- [{csv.name}]({_rel(csv)})")
    md.append("")

    out.write_text("\n".join(md))
    print(f"[p7] report written -> {out}")
    return out


def main():
    build_report()
    print("[p7] reporting complete.")


if __name__ == "__main__":
    main()
