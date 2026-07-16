"""
Global plotnine theme so every figure has a consistent ggplot2 / R look.

Import ``apply_ggplot_theme()`` once at the top of any phase script, or rely on
the module-level side effect (importing this module sets the theme).
"""
from __future__ import annotations

from plotnine import (
    element_line,
    element_rect,
    element_text,
    theme,
    theme_bw,
    theme_set,
)

# A tuned theme_bw that reads like a clean ggplot2 figure.
GGPLOT_THEME = theme_bw(base_size=12) + theme(
    figure_size=(7, 5),
    dpi=150,
    panel_grid_minor=element_line(size=0.25, color="#e8e8e8"),
    panel_grid_major=element_line(size=0.4, color="#dcdcdc"),
    panel_border=element_rect(color="#4d4d4d", size=0.6, fill=None),
    axis_title=element_text(weight="bold"),
    plot_title=element_text(weight="bold", size=14),
    legend_key=element_rect(fill="white", color="white"),
    strip_background=element_rect(fill="#ececec", color="#4d4d4d"),
    strip_text=element_text(weight="bold", size=10),
)

# A curated, colour-blind-friendly palette (Okabe-Ito) for discrete groups.
OKABE_ITO = [
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # green
    "#CC79A7",  # purple/pink
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
]


def apply_ggplot_theme() -> None:
    """Set the global plotnine theme. Call once per script (idempotent)."""
    theme_set(GGPLOT_THEME)


# Side effect on import: make the theme active immediately.
apply_ggplot_theme()
