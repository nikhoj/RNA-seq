"""
Warmup pipeline (airway / GSE52778): Phases 0 -> 3 + 5 + 7.

No time axis, so Phase 4 (time-course) and the time parts of Phase 6 are skipped.
Runs everything end-to-end and writes figures, tables, and a report.

Run:  python scripts/run_warmup.py
"""
import os
import sys
from pathlib import Path

os.environ["RNASEQ_DATASET"] = "airway"
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import (  # noqa: E402
    p0_setup, p1_download, p1_wrangle, p2_qc, p3_de, p5_enrichment, p7_report,
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
    print("\n########## PHASE 5: enrichment ##########")
    p5_enrichment.main()
    print("\n########## PHASE 7: report ##########")
    p7_report.main()
    print("\n=== WARMUP PIPELINE COMPLETE ===")


if __name__ == "__main__":
    main()
