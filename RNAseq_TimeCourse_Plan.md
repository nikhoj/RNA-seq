# RNA-seq Time-Course Analysis & Visualization — Complete Learning Plan

> **Goal:** Learn bulk RNA-seq analysis end-to-end using a time-course (drug/disease-over-time) dataset, with **all visualizations done in Python but styled to look like R / ggplot2**.

---

## Datasets

### Warmup dataset (1–2 days) — `airway` / GSE52778
- Dexamethasone-treated human airway smooth muscle cells (4 cell lines, treated vs untreated, 18 h).
- **No time axis** — but it is the most-documented RNA-seq teaching dataset in existence.
- Installs in seconds via PyDESeq2 (`from pydeseq2.utils import load_example_data`) — zero download headache.
- Purpose: learn the machinery (counts → normalization → DE → volcano) risk-free.

### Main project — GSE212041 (COVID-19 neutrophil longitudinal cohort)
- Longitudinal **human** samples: 306 hospitalized COVID-19+ patients, 78 symptomatic controls, 8 healthy controls.
- **781 samples across multiple time points** — real temporal axis.
- Processed data provided as **raw counts + TPM** in GEO Supplementary files (ready-made matrix, no alignment needed).
- Purpose: the main event — everything time-course-specific.
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212041

> **Cleaner alternative** if you want a tidier *drug* design: any differentiation/treatment series sampled at day 0/2/4/6 (balanced, ~3 replicates/group). Fewer analyses, easier bookkeeping.

---

## The Python "looks-like-R" stack

| Purpose | Package | R equivalent |
|---|---|---|
| Differential expression | **PyDESeq2** | DESeq2 |
| ggplot2-style plotting | **plotnine** | ggplot2 |
| Data container + PCA/UMAP/clustering | **scanpy / anndata** | SummarizedExperiment |
| Enrichment / GSEA | **gseapy** | clusterProfiler / fgsea |
| Heatmaps, fallback plots | **seaborn / matplotlib** | pheatmap / ComplexHeatmap |
| Volcano label placement | **adjustText** | ggrepel |
| Core numerics | pandas, numpy, scipy, statsmodels | base R / dplyr |

Install:
```bash
pip install pydeseq2 plotnine scanpy anndata gseapy seaborn adjusttext \
            pandas numpy scipy statsmodels GEOparse
```

---

## Phase 0 — Environment & setup

- [ ] Create a clean virtual environment (conda or venv).
- [ ] Install the stack above.
- [ ] Verify imports and versions in a notebook.
- [ ] Set a global plotnine theme so every figure has a consistent ggplot look:
  ```python
  from plotnine import theme_bw, theme_set
  theme_set(theme_bw())
  ```

**Deliverable:** a notebook that imports everything and prints versions.

---

## Phase 1 — Data acquisition & wrangling

- [ ] **Warmup:** load `airway` counts + metadata directly from PyDESeq2.
- [ ] **Main:** download GSE212041 supplementary count matrix + sample metadata (via `GEOparse` or manual download from GEO).
- [ ] Build a clean **counts (genes × samples)** table and a **metadata (samples × attributes)** table.
- [ ] Align sample IDs between counts and metadata; check they match exactly.
- [ ] Parse the **time** and **condition** columns into tidy factors.
- [ ] Filter low-count genes (e.g. keep genes with ≥10 counts in ≥N samples).

**Visuals:** library-size bar plots, detected-genes-per-sample, count-distribution densities (log scale).

---

## Phase 2 — Exploratory data analysis (QC)

- [ ] Normalize with DESeq2 median-of-ratios (size factors).
- [ ] Apply variance-stabilizing transform (VST) for visualization.
- [ ] Run PCA on the VST matrix.

**Visuals:**
- PCA scatter colored by **time**, shaped by **condition** (plotnine).
- Sample-to-sample distance heatmap (seaborn clustermap).
- Hierarchical clustering dendrogram.
- Per-sample normalized-count boxplots (before vs after normalization).
- MA plots.

---

## Phase 3 — Differential expression (core)

- [ ] Build the DESeq2 dataset with an appropriate design formula (e.g. `~ condition` or `~ patient + time`).
- [ ] Run pairwise contrasts:
  - each timepoint vs baseline, and/or
  - disease vs control.
- [ ] Extract results tables (log2FC, p-value, padj).
- [ ] Apply shrinkage to log2 fold changes.
- [ ] Define significance thresholds (e.g. padj < 0.05, |log2FC| ≥ 1).

**Visuals:**
- **Volcano plots** (plotnine + adjustText for gene labels).
- MA plots per contrast.
- Top-DE-gene **normalized-count boxplots / beeswarm** (the classic ggplot look).
- Bar summary of up/down DE gene counts per contrast.

---

## Phase 4 — Time-course-specific analysis (the main reason)

This is what snapshot datasets cannot teach.

- [ ] **Time-series-aware DE (LRT):** likelihood-ratio test comparing a full model (`~ condition + time + condition:time`) vs a reduced model (`~ condition + time`) to find genes whose trajectory differs over time.
- [ ] **Temporal clustering:** cluster significant genes by trajectory *shape* (k-means, or soft/fuzzy clustering — the Python equivalent of Mfuzz).
- [ ] **Spline / trajectory modeling:** concepts behind maSigPro & ImpulseDE2 — fit smooth expression curves over time.
- [ ] Characterize each cluster (early responders, late responders, transient, monotonic up/down).

**Visuals:**
- Expression-trajectory **line plots** (mean ± ribbon per group) with plotnine.
- Cluster-averaged trajectory **panels** (facet_wrap).
- Time-ordered **heatmaps** (genes × timepoints), rows clustered.
- Ribbon/confidence-band plots for individual genes of interest.

---

## Phase 5 — Functional / pathway analysis

- [ ] Over-representation analysis (ORA) on DE gene sets (Enrichr via gseapy).
- [ ] GSEA on ranked gene lists.
- [ ] Analyze enrichment **per temporal cluster** to give each trajectory a biological identity.

**Visuals:**
- Enrichment **dot plots** and bar plots.
- GSEA running-enrichment curves.
- Pathway-level heatmaps across time.
- Enrichment-map / network plots (optional).

---

## Phase 6 — Advanced & integrative (optional)

- [ ] Gene co-expression modules (WGCNA-style correlation network).
- [ ] Module–trait / module–time relationship heatmaps.
- [ ] Signature scoring across time (e.g. score a gene set per sample, plot over time).
- [ ] Optional ML: classify outcome (survivor vs non-survivor) from expression; evaluate with ROC.

**Visuals:** module-trait heatmaps, correlation networks, signature-score-over-time line plots, ROC curves.

---

## Phase 7 — Reporting

- [ ] Assemble a clean notebook / report with the strongest figures.
- [ ] Write a narrative of biological findings per phase.
- [ ] Export publication-style figures (consistent theme, fonts, sizes).
- [ ] Summarize: which genes/pathways change over time, and what it means biologically.

---

## Suggested learning arc

1. **Days 1–2:** Phases 0–3 on the `airway` warmup (fast, well-documented, builds confidence).
2. **Days 3–5:** Repeat Phases 0–3 on GSE212041 (real data, messier metadata).
3. **Days 6+:** Phases 4–7 on GSE212041 — the time-course heart of the project.

---

## Key references

- PyDESeq2 docs — https://pydeseq2.readthedocs.io
- plotnine (ggplot2 for Python) — https://plotnine.org
- DESeq2 RNA-seq workflow (concepts) — https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
- gseapy — https://gseapy.readthedocs.io
- GSE212041 (main dataset) — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212041
- GSE52778 (airway warmup) — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778
