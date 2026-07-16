"""
Phase 1a - Data acquisition.

Both datasets are fetched the SAME way, using GEO's NCBI-generated RNA-seq count
matrices (Entrez GeneID x GSM sample) plus the shared human gene-annotation
table (GeneID -> Symbol). Sample metadata is parsed from the series SOFT file
via GEOparse.

  - airway    (GSE52778)  : 8 samples, dex-treated vs untreated (~740 KB)
  - gse212041 (GSE212041) : 781 longitudinal COVID-19 neutrophil samples

Run standalone to just fetch + cache the raw tables:

    RNASEQ_DATASET=airway    python -m src.p1_download
    RNASEQ_DATASET=gse212041 python -m src.p1_download

Note: PyDESeq2 0.5.x no longer bundles the airway dataset (only synthetic data),
so we download the real airway counts from GEO. This also gives real gene
symbols, which Phase 5 enrichment needs.
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402


# --------------------------------------------------------------------------- #
# Download helpers
# --------------------------------------------------------------------------- #
def _download(url: str, dest: Path) -> Path:
    """Download url -> dest unless already present and non-trivial in size."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"[p1] cached {dest.name} ({dest.stat().st_size:,} bytes)")
        return dest
    print(f"[p1] downloading {dest.name} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
        f.write(r.read())
    print(f"[p1] saved {dest.name} ({dest.stat().st_size:,} bytes)")
    return dest


def _raw_dir() -> Path:
    return config.RAW_AIRWAY if config.DATASET == "airway" else config.RAW_GSE212041


def load_annotation() -> pd.DataFrame:
    """GeneID -> Symbol/GeneType annotation, shared across all human series."""
    dest = _raw_dir() / "Human.GRCh38.p13.annot.tsv.gz"
    _download(config.NCBI_ANNOT_URL, dest)
    annot = pd.read_csv(
        dest, sep="\t", compression="gzip",
        usecols=["GeneID", "Symbol", "GeneType", "EnsemblGeneID"],
    ).set_index("GeneID")
    return annot


# --------------------------------------------------------------------------- #
# GSE212041 author-supplied matrix (Ensembl Gene.ID + Symbol column)
# --------------------------------------------------------------------------- #
def load_gse212041_counts() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Download the author count matrix and split it into:
      - counts : genes (Ensembl Gene.ID) x samples (title, e.g. '1_D0')
      - annot  : Gene.ID -> Symbol
    """
    dest = config.RAW_GSE212041 / "GSE212041_Count_Matrix.txt.gz"
    _download(config.GSE212041_COUNTS_URL, dest)
    df = pd.read_csv(dest, sep="\t", compression="gzip")
    df = df.rename(columns={"Gene.ID": "gene"}).set_index("gene")
    annot = df[["Symbol"]].copy()
    counts = df.drop(columns=["Symbol"])
    counts.index.name = "gene"
    print(f"[p1] gse212041 counts {counts.shape} (Ensembl genes x samples)")
    return counts, annot


# --------------------------------------------------------------------------- #
# Counts + metadata
# --------------------------------------------------------------------------- #
def load_ncbi_counts() -> pd.DataFrame:
    """Raw-count matrix (genes x samples) with Entrez GeneID index."""
    acc = config.ACCESSION[config.DATASET]
    dest = _raw_dir() / f"{acc}_raw_counts_NCBI.tsv.gz"
    _download(config.NCBI_COUNTS_URL.format(acc=acc), dest)
    counts = pd.read_csv(dest, sep="\t", index_col=0, compression="gzip")
    counts.index.name = "gene"
    print(f"[p1] counts {counts.shape} (genes x samples)")
    return counts


def load_soft_metadata(index_by: str = "gsm") -> pd.DataFrame:
    """
    Parse per-sample metadata from the GEO series SOFT file.

    index_by="gsm"   -> index is the GSM accession (airway count columns are GSMs)
    index_by="title" -> index is the sample title (GSE212041 columns are titles)
    """
    import GEOparse

    acc = config.ACCESSION[config.DATASET]
    gse = GEOparse.get_GEO(
        geo=acc, destdir=str(_raw_dir()),
        annotate_gpl=False, include_data=False, silent=True,
    )
    rows: dict[str, dict] = {}
    for gsm_name, gsm in gse.gsms.items():
        title = gsm.metadata.get("title", [""])[0]
        attrs: dict[str, object] = {
            "gsm": gsm_name,
            "title": title,
            "source": gsm.metadata.get("source_name_ch1", [""])[0],
        }
        # characteristics_ch1 holds "key: value" pairs -> flatten to columns
        for item in gsm.metadata.get("characteristics_ch1", []):
            if ":" in item:
                key, value = item.split(":", 1)
                col = key.strip().lower().replace(" ", "_").replace("-", "")
                attrs[col] = value.strip()
        key = title if index_by == "title" else gsm_name
        rows[key] = attrs
    meta = pd.DataFrame.from_dict(rows, orient="index")
    meta.index.name = "sample"
    print(f"[p1] metadata {meta.shape}; columns: {list(meta.columns)}")
    return meta


def acquire() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (counts, metadata, annotation) for the active dataset."""
    if config.DATASET == "gse212041":
        counts, annot = load_gse212041_counts()
        meta = load_soft_metadata(index_by="title")
    else:
        counts = load_ncbi_counts()
        meta = load_soft_metadata(index_by="gsm")
        annot = load_annotation()
    # Keep metadata rows in the same order the counts columns appear.
    shared = [s for s in counts.columns if s in meta.index]
    meta = meta.loc[shared]
    return counts, meta, annot


def main() -> None:
    from src import io_utils

    counts, meta, annot = acquire()
    io_utils.save_counts(counts, "counts_raw.parquet")
    io_utils.save_metadata(meta, "metadata_raw.parquet")
    annot.to_parquet(config.processed_path("annotation.parquet"))
    print("[p1] raw acquisition cached.")


if __name__ == "__main__":
    main()
