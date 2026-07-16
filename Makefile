# RNA-seq time-course pipeline
# Uses the `rnaseq` conda env's python. Override with:  make PY=python warmup
PY ?= /Users/abdullah/miniforge3/envs/rnaseq/bin/python

.PHONY: help setup warmup main clean-figs clean-processed

help:
	@echo "make setup          - verify environment (Phase 0)"
	@echo "make warmup         - run the airway warmup pipeline (Phases 0-3,5,7)"
	@echo "make main           - run the full GSE212041 pipeline (Phases 0-7)"
	@echo "make clean-figs      - delete generated figures"
	@echo "make clean-processed - delete cached processed data"

setup:
	RNASEQ_DATASET=airway $(PY) -m src.p0_setup

warmup:
	$(PY) scripts/run_warmup.py

main:
	$(PY) scripts/run_main.py

# --- individual airway phases (handy while iterating) ---------------------
p1:
	RNASEQ_DATASET=airway $(PY) -m src.p1_download && RNASEQ_DATASET=airway $(PY) -m src.p1_wrangle
p2:
	RNASEQ_DATASET=airway $(PY) -m src.p2_qc
p3:
	RNASEQ_DATASET=airway $(PY) -m src.p3_de
p5:
	RNASEQ_DATASET=airway $(PY) -m src.p5_enrichment

clean-figs:
	rm -rf results/figures/*/*.png

clean-processed:
	rm -f data/processed/*.parquet
