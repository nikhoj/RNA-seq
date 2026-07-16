"""
Main pipeline (GSE212041 COVID-19 neutrophil longitudinal cohort): Phases 0 -> 7.

This is the full time-course project. Before the first run, inspect the real
SOFT metadata columns and adjust config.DESIGN["gse212041"] (condition_col,
condition_baseline, time_col) to match. The wrangler auto-parses a numeric time
from the configured time_col.

Run:  python scripts/run_main.py
"""
import os
import sys
from pathlib import Path

os.environ["RNASEQ_DATASET"] = "gse212041"
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import (  # noqa: E402
    p0_setup, p1_download, p1_wrangle, p2_qc, p3_de,
    p4_timecourse, p5_enrichment, p6_advanced, p7_report,
)


def main():
    print("\n########## PHASE 0: setup ##########")
    p0_setup.main()
    print("\n########## PHASE 1a: download ##########")
    p1_download.main()
    print("\n########## PHASE 1b: wrangle + QC ##########")
    p1_wrangle.main()
    print("\n########## PHASE 2: EDA / normalization ##########")
    p2_qc.main()
    print("\n########## PHASE 3: differential expression ##########")
    p3_de.main()
    print("\n########## PHASE 4: time-course ##########")
    p4_timecourse.main()
    print("\n########## PHASE 5: enrichment ##########")
    p5_enrichment.main()
    print("\n########## PHASE 6: advanced ##########")
    p6_advanced.main()
    print("\n########## PHASE 7: report ##########")
    p7_report.main()
    print("\n=== MAIN PIPELINE COMPLETE ===")


if __name__ == "__main__":
    main()
