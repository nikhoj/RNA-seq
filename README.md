# RNA-seq Time-Course Analysis (Python, ggplot-styled)

End-to-end bulk RNA-seq analysis — download → QC → differential expression →
time-course modelling → pathway enrichment → report — implemented entirely in
Python with R/ggplot2-style figures. See
[`RNAseq_TimeCourse_Plan.md`](RNAseq_TimeCourse_Plan.md) for the learning plan
this project implements.

Two datasets share one pipeline via the `RNASEQ_DATASET` switch:

| Dataset | Accession | Role | Time axis |
|---|---|---|---|
| `airway` | GSE52778 | warmup (dex vs untreated, 8 samples) | no |
| `gse212041` | GSE212041 | main COVID-19 neutrophil cohort (781 samples) | yes |

## Setup

```bash
conda create -n rnaseq python=3.11 -y
conda activate rnaseq
pip install -r requirements.txt
```

## Run

```bash
make setup     # Phase 0: verify the environment
make warmup    # airway: Phases 0-3, 5, 7  (fully working, ~2 min)
make main      # GSE212041: Phases 0-7      (see note below)
```

Or run any phase standalone (the `RNASEQ_DATASET` env var selects the dataset):

```bash
RNASEQ_DATASET=airway python -m src.p3_de
```

## Pipeline phases

| Phase | Module | Output |
|---|---|---|
| 0 setup | `src/p0_setup.py` | version check, global ggplot theme |
| 1 download | `src/p1_download.py` | GEO NCBI counts + annotation + SOFT metadata |
| 1 wrangle | `src/p1_wrangle.py` | filtered counts, tidy factors, QC plots |
| 2 EDA | `src/p2_qc.py` | size factors, VST, PCA, distance heatmap |
| 3 DE | `src/p3_de.py` | Wald contrasts, volcano, MA, top-gene boxplots |
| 4 time-course | `src/p4_timecourse.py` | LRT, fuzzy trajectory clusters, time heatmaps |
| 5 enrichment | `src/p5_enrichment.py` | Enrichr ORA + pre-ranked GSEA |
| 6 advanced | `src/p6_advanced.py` | co-expression modules, signatures, outcome ROC |
| 7 report | `src/p7_report.py` | `report/<dataset>_report.md` |

Outputs land in `results/figures/`, `results/tables/`, and `report/`.
Intermediate state is cached in `data/processed/` (parquet), so later phases
never re-download or re-fit.

## Layout

```
config.py              # dataset switch, paths, thresholds, per-dataset design
src/                   # theme + io_utils + one module per phase
scripts/               # run_warmup.py (airway), run_main.py (GSE212041)
data/raw|processed/    # downloads + cached tables  (gitignored)
results/figures|tables # generated outputs          (gitignored)
report/                # generated Markdown report
```

## Running the main dataset (GSE212041)

The warmup is validated end-to-end. For the 781-sample main cohort, before the
first `make main`:

1. Run `RNASEQ_DATASET=gse212041 python -m src.p1_download` to fetch counts +
   metadata.
2. Inspect the SOFT metadata columns it prints, then set
   `config.DESIGN["gse212041"]` — `condition_col`, `condition_baseline`,
   `time_col` — to match the real column names.
3. `make main`.

## Notes on the Python "looks-like-R" stack

- **PyDESeq2** for DESeq2 median-of-ratios normalization, VST, and Wald
  contrasts. (v0.5.x no longer bundles airway data, so we download the real
  counts from GEO — this also gives gene symbols for enrichment.)
- **LRT for time-course** is implemented via `statsmodels` (full vs reduced
  model on the VST matrix), since PyDESeq2 0.5.x has no built-in LRT.
- **Fuzzy trajectory clustering** uses `scikit-fuzzy` (the Python Mfuzz), with a
  k-means fallback.
- **plotnine** for ggplot2 figures, **adjustText** for volcano labels,
  **gseapy** for Enrichr/GSEA, **seaborn** for clustermaps.
