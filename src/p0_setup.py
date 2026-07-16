"""
Phase 0 - Environment & setup.

Verifies every package in the stack imports, prints versions, sets the global
plotnine theme, and creates the output directory tree.

Run:  python -m src.p0_setup
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402

PACKAGES = [
    "pydeseq2",
    "plotnine",
    "scanpy",
    "anndata",
    "gseapy",
    "seaborn",
    "adjustText",
    "pandas",
    "numpy",
    "scipy",
    "statsmodels",
    "GEOparse",
    "sklearn",
    "matplotlib",
]


def main() -> None:
    print("=" * 60)
    print("Phase 0 - environment check")
    print("=" * 60)
    print(f"Python: {sys.version.split()[0]}")
    ok, missing = [], []
    for name in PACKAGES:
        try:
            mod = importlib.import_module(name)
            ver = getattr(mod, "__version__", "?")
            ok.append(name)
            print(f"  [OK]   {name:<14} {ver}")
        except Exception as exc:  # noqa: BLE001
            missing.append(name)
            print(f"  [MISS] {name:<14} ({exc})")

    # Set the global ggplot theme (side effect of importing the module).
    from src import theme  # noqa: F401

    config.ensure_dirs()
    print("-" * 60)
    print(f"Imported {len(ok)}/{len(PACKAGES)} packages; theme set; dirs ready.")
    if missing:
        print(f"MISSING: {', '.join(missing)} -- install before proceeding.")
        sys.exit(1)
    print("Phase 0 complete.")


if __name__ == "__main__":
    main()
