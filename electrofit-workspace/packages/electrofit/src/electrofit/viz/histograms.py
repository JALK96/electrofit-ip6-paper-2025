"""Histogram plotting utilities for charge ensemble diagnostics.

Refactored from legacy step6 script (``6_extract_charges.py``) and
``analysis/plot_charges.py`` so that optional diagnostic visualisations
are available in the new pipeline without retaining monolithic code.

Usage patterns (all optional):
 - Per‑atom histograms before any filtering (``hist.pdf``)
 - After outlier removal (``hist_no_outlier.pdf``)
 - After group/symmetry averaging (``hist_adjusted.pdf``)

Public functions here are intentionally low level; orchestration / naming
decisions are performed inside the domain layer (``average_charges.py``).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import matplotlib.pyplot as plt  # Hard dependency – removal of .skip sentinel logic

__all__ = [
    "plot_atom_histograms",
    "HistogramSpec",
]


@dataclass
class HistogramSpec:
    """Specification for a single histogram panel.

    Attributes
    ----------
    column : str
        Atom (column) name.
    data_before : Sequence[float]
        Raw (pre-filter) charges.
    data_after : Sequence[float] | None
        Optional second dataset (e.g. post-filter, post-average) for overlay.
    mean_line : float | None
        Value at which to draw vertical mean line (typically adjusted/group mean).
    color : str
        Base color for the column / group.
    overlay_color : str | None
        Colour for overlay dataset (falls back to base color with lower alpha).
    combine_group : bool
        If True and *group_atoms* provided externally, panels may use group-combined
        distributions (implemented upstream by passing pre-combined data).
    """

    column: str
    data_before: Sequence[float]
    data_after: Sequence[float] | None = None
    mean_line: float | None = None
    color: str = "darkblue"
    overlay_color: str | None = None
    combine_group: bool = False


def _auto_bins(a: Sequence[float], b: Sequence[float] | None, bins: int):
    data = list(a)
    if b:
        data.extend(b)
    if not data:
        return np.linspace(-1, 1, 5)  # fallback
    return np.histogram_bin_edges(data, bins=bins)


def plot_atom_histograms(
    specs: Sequence[HistogramSpec],
    outfile: str | Path,
    title: str,
    bins: int = 20,
    ncols: int = 4,
) -> None:
    """Render a grid of histograms based on provided :class:`HistogramSpec`s.

    Parameters
    ----------
    specs : list[HistogramSpec]
        One spec per panel. Order preserved.
    outfile : str | Path
        Output PDF path.
    title : str
        Figure suptitle.
    bins : int
        Target number of bins (used to derive common bin edges per panel).
    ncols : int
        Maximum columns per row (layout adapts automatically).
    """
    if not specs:
        return
    if not specs:
        return
    # Hard import already executed; any ImportError surfaces earlier.
    n = len(specs)
    ncols = min(ncols, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    if not isinstance(axes, np.ndarray):  # single panel corner case
        axes = np.array([axes])
    axes = axes.flatten()
    for i, spec in enumerate(specs):
        ax = axes[i]
        data_orig = [v for v in spec.data_before if v is not None]
        data_after = (
            [v for v in spec.data_after if v is not None]
            if spec.data_after is not None
            else None
        )
        edges = _auto_bins(data_orig, data_after, bins)
        # Original / before
        ax.hist(
            data_orig,
            bins=edges,
            color=spec.color,
            alpha=0.5 if data_after is not None else 0.9,
            edgecolor="black",
            label="Original" if data_after is not None else "Distribution",
        )
        if data_after is not None:
            ax.hist(
                data_after,
                bins=edges,
                color=spec.overlay_color or spec.color,
                alpha=0.9,
                edgecolor="black",
                label="Adjusted",
            )
        if spec.mean_line is not None:
            ax.axvline(
                spec.mean_line,
                color="red" if data_after is not None else "black",
                linestyle="dashed",
                linewidth=1.5,
                label=f"{spec.mean_line:.2f}",
            )
        ax.set_title(spec.column)
        ax.set_xlabel("Charge")
        ax.set_ylabel("Frequency")
        ax.legend(fontsize=8)
    # Hide any unused axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile)
    plt.close(fig)
