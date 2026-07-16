"""
Central configuration for the RNA-seq time-course project.

Everything path- and threshold-related lives here so the phase scripts stay
clean. The ``DATASET`` switch lets the *same* pipeline run on the airway warmup
and on the GSE212041 main project.
"""
from __future__ import annotations

import os
from pathlib import Path

# --------------------------------------------------------------------------- #
# Dataset switch
# --------------------------------------------------------------------------- #
# "airway"    -> warmup (Phases 0-3), loaded straight from PyDESeq2, no download
# "gse212041" -> main COVID-19 neutrophil longitudinal cohort (Phases 0-7)
DATASET = os.environ.get("RNASEQ_DATASET", "airway")

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"
REPORT = ROOT / "report"

# Per-dataset raw download folders
RAW_AIRWAY = DATA_RAW / "airway"
RAW_GSE212041 = DATA_RAW / "gse212041"

# GEO accessions
GSE_WARMUP = "GSE52778"    # airway (dexamethasone) warmup
GSE_MAIN = "GSE212041"     # COVID-19 neutrophil longitudinal cohort

# Per-dataset accession lookup
ACCESSION = {"airway": GSE_WARMUP, "gse212041": GSE_MAIN}

# NCBI-generated RNA-seq count download endpoints (used for the airway warmup)
NCBI_COUNTS_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?type=rnaseq_counts&acc={acc}&format=file&file={acc}_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)
NCBI_ANNOT_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
)

# GSE212041 ships author-supplied supplementary matrices (no NCBI-generated set).
# Row index is Ensembl Gene.ID with a Symbol column; sample columns match the
# GEO sample *title* (e.g. "1_D0"), not the GSM accession.
GSE212041_COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212041/suppl/"
    "GSE212041_Neutrophil_RNAseq_Count_Matrix.txt.gz"
)

# --------------------------------------------------------------------------- #
# Analysis thresholds
# --------------------------------------------------------------------------- #
# Low-count gene filter: keep genes with >= MIN_COUNT reads in >= MIN_SAMPLES samples
MIN_COUNT = 10
MIN_SAMPLES = 3

# Differential-expression significance
PADJ_THRESHOLD = 0.05
LFC_THRESHOLD = 1.0          # |log2 fold change|

# Temporal clustering
N_TEMPORAL_CLUSTERS = 6

# Reproducibility
RANDOM_STATE = 42

# --------------------------------------------------------------------------- #
# Per-dataset design specification
# --------------------------------------------------------------------------- #
# Tells the wrangling step which raw metadata columns map to the analysis
# factors, the baseline level, and (for time-course) the time column.
DESIGN = {
    "airway": {
        "condition_col": "treatment",
        "condition_baseline": "Untreated",
        "condition_keep": ["Untreated", "Dexamethasone"],  # canonical airway subset
        "block_col": "cell_line",     # paired design: ~ cell_line + condition
        "time_col": None,             # no time axis in the warmup
        "design_factors": ["cell_line", "condition"],
    },
    "gse212041": {
        # Real GSE212041 SOFT columns (flattened from characteristics_ch1).
        # Main contrast: COVID+ vs COVID- symptomatic across days 0/3/7.
        "condition_col": "patient_category",
        "condition_baseline": "covid-_symptomatic",
        # keep the two symptomatic arms with a real time course; 'healthy'
        # (n=8, single 'H' timepoint) is excluded from the DE/time-course model
        "condition_keep": ["COVID+", "COVID- symptomatic"],
        "block_col": None,
        "time_col": "time_point",     # D0/D3/D7 -> 0/3/7; DE/H -> NaN (dropped from Phase 4)
        "design_factors": ["condition"],
        "outcome_col": "acuity.max",  # severity 1-5, used by Phase 6 ML
        # Time-course test. Controls are D0-only here, so a condition:time
        # interaction is unidentifiable; the meaningful question is which genes
        # change over the disease course *within* COVID+.
        "timecourse_test": "time_main",       # "interaction" | "time_main"
        "timecourse_condition": "covid+",     # restrict the trajectory analysis
        # Binarise acuity for the Phase 6 outcome classifier:
        #   severe = acuity 1-2 (intubated/deceased), mild = acuity 4-5.
        "outcome_positive": ["1", "2"],
        "outcome_negative": ["4", "5"],
    },
}


def design() -> dict:
    """Design spec for the active dataset."""
    return DESIGN[DATASET]

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def fig_dir(phase: str) -> Path:
    """Return (and create) the figure directory for a given phase folder name."""
    d = FIGURES / phase
    d.mkdir(parents=True, exist_ok=True)
    return d


def ensure_dirs() -> None:
    """Create every output directory if missing (idempotent)."""
    for d in (
        DATA_RAW, DATA_PROCESSED, RESULTS, TABLES, FIGURES, REPORT,
        RAW_AIRWAY, RAW_GSE212041,
    ):
        d.mkdir(parents=True, exist_ok=True)


def processed_path(name: str) -> Path:
    """Path to a processed artifact, namespaced by the active dataset."""
    DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
    return DATA_PROCESSED / f"{DATASET}_{name}"


if __name__ == "__main__":
    ensure_dirs()
    print(f"Active dataset : {DATASET}")
    print(f"Project root   : {ROOT}")
    print(f"Processed dir  : {DATA_PROCESSED}")
    print(f"Figures dir    : {FIGURES}")
