from __future__ import annotations
import logging
import pathlib
from typing import Dict, Tuple, Literal

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array         #  (https://userguide.mdanalysis.org/examples/analysis/distances_and_contacts/distances_between_selections.html)
from MDAnalysis.analysis.rdf import InterRDF

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import re
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from typing import Optional
from mpl_toolkits.mplot3d import Axes3D
sns.set_context("talk")

RCUTOFF_NM         = 0.3       # coordination threshold  (change if desired)

# Mapping: phosphate → peripheral O atom names (no ester O)
PHOS_OXYGENS: Dict[str, Tuple[str, ...]] = {
    "P" : ("O6",  "O7",  "O8"),
    "P1": ("O9",  "O10", "O11"),
    "P2": ("O12", "O13", "O14"),
    "P3": ("O15", "O16", "O17"),
    "P4": ("O18", "O19", "O20"),
    "P5": ("O21", "O22", "O23"),
}
PHOS_LABELS = tuple(PHOS_OXYGENS.keys())   # ('P', 'P1', …)
PHOS_LABELS_CORRECTED = ('P1', 'P2', 'P3', 'P4', 'P5', 'P6')  # re-order to match the paper


# ╔═══════════════════════════════════════════════════════════════╗
# ┃                         PLOTTING                              ┃
# ╚═══════════════════════════════════════════════════════════════╝

def plot_counts_subplots(
        counts_ts : np.ndarray,          # shape (6, n_frames)
        timestep_ps: float,
        out_png   : pathlib.Path,
        title     : str,
        *,
        ion_label: str = "Na⁺",
) -> None:
    """Re-implementation of your 2×3 layout but plotting counts."""
    sns.set_context("talk")
    n_phos, n_frames = counts_ts.shape
    time_ns = np.arange(n_frames) * timestep_ps / 1000.0

    fig, axes = plt.subplots(2, 3, figsize=(8, 6), sharex=True, sharey=True)
    axes = axes.flatten()

    for idx, (ax, label) in enumerate(zip(axes, PHOS_LABELS)):
        y = counts_ts[idx]
        ax.plot(time_ns, y, lw=0.7, color="black")
        mean_y = y.mean()
        ax.axhline(mean_y, ls="--", lw=1.4, color="red")
        ax.text(0.95, 0.9, f"⟨count⟩ = {mean_y:.1f}",
                ha="right", va="top", transform=ax.transAxes, color="red", fontsize=11)
        ax.set_title(f"{label} – {ion_label}")

        if idx % 3 == 0:
            ax.set_ylabel(f"No. bound {ion_label}")
        if idx >= 3:
            ax.set_xlabel("Time / ns")

    fig.suptitle(title, y=1.02, fontsize=16)
    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logging.info("Plot written to %s", out_png)


def plot_coordination_boxplot(
        counts_ts: np.ndarray,          # shape (6, n_frames)
        out_png: pathlib.Path,
        *,
        showfliers: bool = True,
        binary_code: Optional[str] = None,
        ylabel: str = "Na⁺ count",
        ion_label: str = "Na⁺",
) -> None:
    """
    Boxplot of Na⁺ coordination counts per phosphate for one microstate.

    Styling and colour cues are lifted from the summarize-nap boxplots:
    boxes corresponding to a "1" bit in the microstate code are shaded blue.
    """
    if counts_ts.ndim != 2 or counts_ts.shape[0] != len(PHOS_LABELS):
        raise ValueError("counts_ts must have shape (6, n_frames)")

    # infer binary code from output path if not provided explicitly
    if binary_code is None:
        try:
            parent_name = out_png.parents[2].name  # e.g. IP_101101 or IP_101101_100
            match = re.search(r"[01]{6}", parent_name)
            binary_code = match.group(0) if match else "000000"
        except Exception:
            binary_code = "000000"
    if len(binary_code) != 6 or any(c not in "01" for c in binary_code):
        raise ValueError(f"Unexpected binary code '{binary_code}' in output path")
    colours = ["blue" if bit == "1" else "black" for bit in binary_code]

    per_phos = [np.asarray(counts_ts[i], dtype=float) for i in range(len(PHOS_LABELS))]

    metrics = coordination_boxplot_metrics(counts_ts)
    ymin = metrics["ymin"]
    ymax = metrics["ymax"]

    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    bp = ax.boxplot(
        per_phos,
        positions=range(1, 7),
        widths=0.6,
        showfliers=showfliers,
        patch_artist=True,
        medianprops={"visible": False},
        showmeans=True,
        meanprops={
            "marker": "D",
            "markeredgecolor": "black",
            "markerfacecolor": "black",
            "markersize": 6,
        },
    )
    # edge and face colouring per phosphate bit
    for idx, box in enumerate(bp["boxes"]):
        box.set_edgecolor("black")
        if colours[idx] == "blue":
            box.set_facecolor(to_rgba("blue", alpha=0.3))
        else:
            box.set_facecolor("white")
    for element in ("whiskers", "caps"):
        for obj in bp[element]:
            obj.set_color("black")

    ax.set_xticks(range(1, 7))
    ax.set_xticklabels(PHOS_LABELS_CORRECTED, fontsize=12)
    ax.set_ylabel(ylabel.replace("Na⁺", ion_label), fontsize=12)
    ax.set_ylim(ymin, ymax)
    ax.tick_params(axis="y", labelsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(binary_code, fontsize=14)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logging.info("Coordination boxplot written to %s", out_png)


def coordination_boxplot_metrics(counts_ts: np.ndarray) -> Dict[str, object]:
    """Compute the same per-phosphate boxplot stats used for plotting.

    Returns a dict with:
      - per_phosphate: {label: {mean,std,median,q1,q3,iqr,whisker_low,whisker_high,n_outliers}}
      - ymin/ymax: global y-limits based on whiskers (with padding)
    """
    if counts_ts.ndim != 2 or counts_ts.shape[0] != len(PHOS_LABELS):
        raise ValueError("counts_ts must have shape (6, n_frames)")

    per_phosphate: Dict[str, Dict[str, float | int]] = {}
    whisk_lows: list[float] = []
    whisk_highs: list[float] = []

    for idx, label in enumerate(PHOS_LABELS):
        arr = np.asarray(counts_ts[idx], dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            per_phosphate[label] = {
                "mean": float("nan"),
                "std": float("nan"),
                "median": float("nan"),
                "q1": float("nan"),
                "q3": float("nan"),
                "iqr": float("nan"),
                "whisker_low": float("nan"),
                "whisker_high": float("nan"),
                "n_outliers": 0,
            }
            continue

        q1, median, q3 = np.percentile(arr, (25, 50, 75))
        iqr = q3 - q1
        low_bound = q1 - 1.5 * iqr
        high_bound = q3 + 1.5 * iqr

        lw_cand = arr[arr >= low_bound]
        hw_cand = arr[arr <= high_bound]
        whisk_low = float(lw_cand.min()) if lw_cand.size else float(arr.min())
        whisk_high = float(hw_cand.max()) if hw_cand.size else float(arr.max())
        n_outliers = int(np.count_nonzero((arr < low_bound) | (arr > high_bound)))

        whisk_lows.append(whisk_low)
        whisk_highs.append(whisk_high)

        per_phosphate[label] = {
            "mean": float(np.mean(arr)),
            "std": float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0,
            "median": float(median),
            "q1": float(q1),
            "q3": float(q3),
            "iqr": float(iqr),
            "whisker_low": float(whisk_low),
            "whisker_high": float(whisk_high),
            "n_outliers": n_outliers,
        }

    if whisk_lows and whisk_highs:
        ymin, ymax = float(min(whisk_lows)), float(max(whisk_highs))
        pad = 0.05 * (ymax - ymin if ymax > ymin else 1.0)
        ymin -= pad
        ymax += pad
    else:
        ymin, ymax = 0.0, 1.0

    return {"per_phosphate": per_phosphate, "ymin": ymin, "ymax": ymax}


# ┃                 RDF  of Na⁺ vs. peripheral O                  ┃
# ╚═══════════════════════════════════════════════════════════════╝

# Utility: first‑shell boundary from an RDF curve
def first_shell_end(
    r: np.ndarray,
    g: np.ndarray,
    *,
    r_limit_nm: float = 0.8,
    peak_window_nm: Tuple[float, float] = (0.12, 0.45),
    min_window_nm: Tuple[float, float] = (0.18, 0.65),
) -> float:
    """
    Return the radius r (nm) at which the first coordination shell ends,
    i.e. the first *local* minimum of g(r) that follows the first *local*
    maximum (peak).

    This is a heuristic detector intended for noisy RDF curves. In rare
    cases (especially with sparse sampling) the RDF may show small
    oscillations at very short distances that can be mistaken for the first
    peak. To avoid returning unphysical cutoffs near r≈0, the peak search is
    restricted to a physically plausible window (default 0.12–0.45 nm) and
    the minimum is restricted to a plausible post-peak window (default
    0.18–0.65 nm, capped by ``r_limit_nm``). If no clear peak→minimum pair
    can be detected, NaN is returned so the caller can fall back to a fixed
    cutoff.

    Notes
    -----
    *Local* extrema are detected by finite‑difference sign changes:
      peak   : g[i-1] < g[i] >= g[i+1]
      minimum: g[i-1] > g[i] <= g[i+1]

    The default windows bracket the Na–O first-shell peak and its following
    minimum in all micro-states we observed.
    """
    if r.ndim != 1 or g.ndim != 1:
        raise ValueError("r and g must be one‑dimensional arrays")
    if len(r) != len(g):
        raise ValueError("r and g must have the same length")

    r = np.asarray(r, dtype=float)
    g = np.asarray(g, dtype=float)
    if r.size < 3:
        return np.nan
    if not np.all(np.isfinite(r)) or not np.any(np.isfinite(g)):
        return np.nan

    # Ensure r is monotonically increasing (InterRDF bins are, but keep defensive)
    if np.any(np.diff(r) <= 0):
        order = np.argsort(r)
        r = r[order]
        g = g[order]

    peak_lo, peak_hi = peak_window_nm
    peak_mask = (r >= peak_lo) & (r <= min(peak_hi, r_limit_nm))
    peak_candidates = np.where(peak_mask)[0]
    if peak_candidates.size == 0:
        return np.nan

    # Use the strongest peak in the physically plausible window.
    peak_idx = int(peak_candidates[np.nanargmax(g[peak_candidates])])
    if not np.isfinite(g[peak_idx]) or g[peak_idx] <= 0:
        return np.nan

    # Now find the first local minimum after the selected peak.
    N = len(r)
    min_lo, min_hi = min_window_nm
    min_hi = min(min_hi, r_limit_nm)

    # Minimum search interval: after the peak, within the min window.
    start = peak_idx + 1
    # find first index where r >= min_lo, but keep it after peak
    if r[start] < min_lo:
        start = int(np.searchsorted(r, min_lo, side="left"))
        start = max(start, peak_idx + 1)
    end = int(np.searchsorted(r, min_hi, side="right")) - 1
    end = min(end, N - 2)
    if start >= end:
        return np.nan

    min_idx: int | None = None
    for j in range(start, end + 1):
        if g[j] < g[j - 1] and g[j] <= g[j + 1]:
            min_idx = j
            break

    # Fallback: if there is no clear local minimum, pick the global minimum
    # in the interval (still gives a sensible cutoff for smooth curves).
    if min_idx is None:
        segment = g[start : end + 1]
        if not np.any(np.isfinite(segment)):
            return np.nan
        min_idx = int(start + np.nanargmin(segment))

    shell_end = float(r[min_idx])
    # Sanity checks: must be after the peak and within the configured windows.
    if not (np.isfinite(shell_end) and (shell_end > float(r[peak_idx]))):
        return np.nan
    if shell_end < min_lo or shell_end > min_hi:
        return np.nan
    # Require a drop after the peak (avoid returning a flat/noisy plateau).
    if np.isfinite(g[min_idx]) and g[min_idx] >= 0.95 * g[peak_idx]:
        return np.nan

    return shell_end

def plot_rdf_periphO_Na(
        u: mda.Universe,
        out_png: pathlib.Path,
        r_max: float = 1.2,      # nm
        nbins: int   = 240,      # → dr ≈ 0.005 nm
        *,
        y_max_global: Optional[float] = None,
        hide_y_label: bool = False,
        hide_y_ticklabels: bool = False,
        show_shell_cutoff: bool = False,
        cutoff_mode: Literal["fixed", "first-shell"] = "fixed",
        r_cut_nm: float = RCUTOFF_NM,
        rdf_results_override: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]] = None,
        ion_name: str = "NA",
        ion_label: Optional[str] = None,
) -> None:
    """Plot Na⁺–peripheral‑oxygen RDFs for **one** microstate.

    Six rows (P1–P6), shared x‑axis and *global* y‑axis that must be supplied
    via *y_max_global* so **all** microstates use exactly the same scale.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Trajectory of the microstate to analyse.
    out_png : pathlib.Path
        Output filename.
    r_max : float, optional
        Upper bound of the RDF in **nm**.
    nbins : int, optional
        Number of histogram bins.
    y_max_global : float, keyword‑only, required
        *Upper* limit of the y‑axis for **all** figures.  Compute this once
        across the full data set and pass it to every call.  If *None*, an
        exception is raised so the user does not forget.
    hide_y_label : bool, optional
        When `True`, omit the phosphate label text from the y-axis.
    hide_y_ticklabels : bool, optional
        When `True`, suppress y tick labels (tick marks remain visible).
    show_shell_cutoff : bool, optional
        Toggle drawing the green first-shell cutoff guide.
    cutoff_mode : {"fixed", "first-shell"}, optional
        Controls both the vertical cutoff guide and the integration limit used
        to compute the coordination number displayed in each subplot.
        - ``fixed`` uses ``r_cut_nm``.
        - ``first-shell`` uses the first minimum of g(r) after the first peak
          (per phosphate); falls back to ``r_cut_nm`` if no minimum is detected.
    r_cut_nm : float, optional
        Fixed cutoff radius in nm (used when ``cutoff_mode="fixed"`` and as a
        fallback when the first-shell end cannot be detected).
    rdf_results_override : dict | None, optional
        Precomputed RDF arrays in the form ``{label: (r_nm, g_raw)}``. When
        provided, the expensive InterRDF computation is skipped and these
        arrays are used directly (``g_raw`` should be the raw RDF, i.e. before
        multiplying by 3 to obtain per-phosphate curves).
    """

    if y_max_global is None:
        raise ValueError("plot_rdf_periphO_Na requires *y_max_global* so all figures share the same scale.")
    ion_label = ion_label or f"{ion_name.upper()}⁺"

    font_rc = {
        "font.family": "serif",
        "font.serif": ["Nimbus Roman", "Nimbus Roman No9 L", "Times New Roman", "Times"],
    }

    with plt.rc_context(font_rc):
        # ---- compute RDFs ---------------------------------------------------
        ion_sel = f"name {ion_name}"
        na_ag = u.select_atoms(ion_sel)
        periph_dict = {
            label: u.select_atoms("resname I* and name " + " ".join(PHOS_OXYGENS[label]))
            for label in PHOS_LABELS
        }

        rdf_results: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        if rdf_results_override is None:
            for label, oxy_ag in periph_dict.items():
                rdf = InterRDF(oxy_ag, na_ag, range=(0, r_max * 10), nbins=nbins)
                rdf.run()
                rdf_results[label] = (rdf.bins / 10.0, rdf.rdf)  # convert Å → nm
        else:
            missing = [lab for lab in PHOS_LABELS if lab not in rdf_results_override]
            if missing:
                raise ValueError(f"rdf_results_override missing labels: {', '.join(missing)}")
            for label in PHOS_LABELS:
                r_vals, g_vals = rdf_results_override[label]
                rdf_results[label] = (np.asarray(r_vals, dtype=float),
                                      np.asarray(g_vals, dtype=float))

        # ---- figure & axes --------------------------------------------------
        sns.set_context("talk")
        n = len(PHOS_LABELS)
        fig, axes = plt.subplots(n, 1, figsize=(6, 1.2 * n), sharex=True, sharey=True)
        #fig, axes = plt.subplots(n, 1, figsize=(4.5, 1.2 * n), sharex=True, sharey=True)

        # ---- binary code → colours -----------------------------------------
        try:
            parent_name = out_png.parents[2].name
            match = re.search(r"[01]{6}", parent_name)
            binary_code = match.group(0) if match else "000000"
        except Exception:
            binary_code = "000000"
        if len(binary_code) != 6 or any(c not in "01" for c in binary_code):
            raise ValueError(f"Unexpected binary code '{binary_code}' in output path")
        colours = ["blue" if bit == "1" else "black" for bit in binary_code]

        # ---- constants for inset -------------------------------------------
        X_LOWER, X_UPPER = RCUTOFF_NM, 1.2     # nm

        # ---- loop over P1–P6 ------------------------------------------------
        shell_ends_nm: list[float] = []
        for ax, (label, col) in zip(axes, zip(PHOS_LABELS, colours)):
            r, g_raw = rdf_results[label]
            g_per_phos = g_raw * 3                    # 3 O atoms → per phosphate

            # --- determine first‑shell boundary -----------------------------
            shell_end = first_shell_end(r, g_per_phos)
            have_shell_end = not np.isnan(shell_end)
            if have_shell_end:
                shell_ends_nm.append(float(shell_end))

            cutoff_nm = r_cut_nm
            if cutoff_mode == "first-shell":
                if have_shell_end:
                    # extra guard against pathological detections
                    if shell_end < 0.15:
                        logging.info(
                            "First-shell cutoff for %s was %.3f nm (too small); falling back to fixed cutoff %.3f nm",
                            label,
                            shell_end,
                            r_cut_nm,
                        )
                    else:
                        cutoff_nm = shell_end
                else:
                    logging.info(
                        "No first-shell minimum detected for %s; falling back to fixed cutoff %.3f nm",
                        label,
                        r_cut_nm,
                    )

            # draw guides ------------------------------------------------------
            if cutoff_mode == "first-shell" and have_shell_end:
                # integration cutoff at the first-shell end
                ax.axvline(shell_end, color="green", ls=":", lw=1.2)
                ax.text(
                    shell_end,
                    0.95 * y_max_global,
                    f"{shell_end:.2f} nm",
                    rotation=90,
                    va="top",
                    ha="right",
                    color="green",
                    fontsize=8,
                )
                logging.info("First-shell end for %s: %.3f nm", label, shell_end)
            elif show_shell_cutoff and have_shell_end:
                # visual guide only (fixed integration cutoff)
                ax.axvline(shell_end, color="green", ls=":", lw=1.0)
                ax.text(
                    shell_end,
                    0.95 * y_max_global,
                    f"{shell_end:.2f} nm",
                    rotation=90,
                    va="top",
                    ha="right",
                    color="green",
                    fontsize=8,
                )
                logging.info("First-shell end for %s: %.3f nm", label, shell_end)

            if col == "blue":
                ax.fill_between(r, g_per_phos, 0.0, color=col, alpha=0.3)
            ax.plot(r, g_per_phos, lw=1.2, color="black")
            if cutoff_mode == "fixed" or not have_shell_end:
                ax.axvline(r_cut_nm, ymax=0.8, color="black", ls="--", lw=1.0)

            # global y‑scale --------------------------------------------------
            ax.set_ylim(0.0, y_max_global)

            # y‑axis label ----------------------------------------------------
            if not hide_y_label:
                ax.set_ylabel(PHOS_LABELS_CORRECTED[PHOS_LABELS.index(label)],
                              rotation=0, labelpad=40, va="center")
            else:
                ax.set_ylabel("")
            ax.yaxis.set_major_locator(plt.MaxNLocator(2))
            if hide_y_ticklabels:
                ax.set_yticklabels([])
                ax.tick_params(labelleft=False)
            if ax is not axes[-1]:
                ax.tick_params(labelbottom=False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            # coordination number ---------------------------------------------
            mask_cn = r <= cutoff_nm
            vol_nm3 = np.prod(u.dimensions[:3]) / 1000.0  # Å³ → nm³
            rho_ion  = len(na_ag) / vol_nm3
            logging.info(f"Volume: {vol_nm3:.2f} nm³, Density of {ion_label}: {rho_ion:.2f} nm⁻³")
            coord_num = 4.0 * np.pi * rho_ion * np.trapezoid(g_per_phos[mask_cn] * r[mask_cn]**2, x=r[mask_cn])
            ax.text(0.03, 0.8, f"{coord_num:.2f}", transform=ax.transAxes,
                    ha="left", va="top", fontsize=20, color=col, alpha=0.7)

            # --------  OPTIONAL INSET for "blue" (bit == 1), "black" (bit ==0) ---------
            if col == "blue" or col == "black":
                # restrict data to the desired x-window
                mask_inset = (r >= X_LOWER) & (r <= X_UPPER)
                r_inset    = r[mask_inset]
                g_inset    = g_per_phos[mask_inset]

                # create inset: upper-right corner (loc='upper right')
                axins = inset_axes(ax,
                                   width = 1.5,   # 35 % of parent
                                   height = 0.5,
                                   loc = "upper right",
                                   borderpad = 1.0)

                # plot and format
                if col == "blue":
                    axins.fill_between(r_inset, g_inset, 0,
                                   color=col, alpha=0.4)
                axins.plot(r_inset, g_inset, lw=1.0, color="black")

                # set limits so curve fills inset tightly
                axins.set_xlim(X_LOWER, X_UPPER)
                y_min, y_max = g_inset.min(), g_inset.max()
                x_min, x_max = r_inset.min(), r_inset.max()
                axins.set_ylim(y_min*0.98, y_max*1.02)

                # --- show only the top (maximum) y tick -----------------------
                axins.set_yticks([y_max])                 # one tick at ymax
                axins.set_yticklabels([f"{y_max:.0f}"], fontsize=18)   # optional: format label
                axins.tick_params(axis="y",               # keep the tick mark
                                which="both",
                                direction="out",
                                left=False,
                                right=False,
                                labelright=False)
                if ax is not axes[0]:
                    axins.tick_params(axis="x",               # hide x ticks/labels
                                    bottom=False,
                                    labelbottom=False)
                else:
                    axins.set_xticks([x_min, x_max])                 # one tick at ymax
                    axins.set_xticklabels([f"{x_min:.2f}", f"{x_max:.2f}"], fontsize=18)   # optional: format label
                    axins.tick_params(axis="x",               # show x ticks/labels
                                    top=False,
                                    labeltop=True,
                                    bottom=False,
                                    labelbottom=False)

            #ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))  # scientific notation for large y-values

        # ---- common labels & title -----------------------------------------
        axes[-1].set_xlabel("r / nm", fontsize=18)
        fig.suptitle(binary_code, fontsize=26)

        plt.tight_layout()
        fig.subplots_adjust(hspace=0.1)
        fig.savefig(out_png, dpi=300)
        plt.close(fig)
        logging.info("RDF figure written to %s", out_png)
        if shell_ends_nm:
            mean_nm = float(np.mean(shell_ends_nm))
            std_nm = float(np.std(shell_ends_nm, ddof=1)) if len(shell_ends_nm) > 1 else 0.0
            logging.info(
                "First-shell end summary (n=%d): mean=%.3f nm, std=%.3f nm",
                len(shell_ends_nm),
                mean_nm,
                std_nm,
            )


# ┃    2-D PROJECTION  PERPENDICULAR  TO  A  REFERENCE PLANE      ┃
# ╚═══════════════════════════════════════════════════════════════╝

def plot_frame_network_plane(
        u: mda.Universe,
        frame: int,
        out_png: pathlib.Path,
        reference_triplet: Tuple[str, str, str] = ("P", "P2", "P3"),
        r_cut_nm: float = RCUTOFF_NM,
        ion_name: str = "NA",
        ion_label: Optional[str] = None,
) -> None:
    """
    Draw a 2-D network diagram of ions and phosphate centres by
    projecting onto the plane defined by *any* three IP6 atoms.

    Parameters
    ----------
    u : MDAnalysis.Universe
        MD trajectory already loaded.
    frame : int
        Index of the frame to visualise.
    out_png : pathlib.Path
        Output filename (.png).
    reference_triplet : (str, str, str)
        Names of *three* distinct IP6 atoms whose plane you wish to view
        (e.g. ("P","O9","P4") or any others in the same residue).
    r_cut_nm : float
        Only draw a dashed line when ion–peripheral‐O distance < r_cut_nm,
        or `None` to draw all links.
    """
    ion_label = ion_label or f"{ion_name.upper()}⁺"
    # 1. frame
    n_frames = len(u.trajectory)
    if not (0 <= frame < n_frames):
        raise IndexError(f"Frame {frame} out of range (0…{n_frames-1})")
    u.trajectory[frame]

    # 2. selections
    na_atoms = u.select_atoms(f"name {ion_name}")
    if len(na_atoms) == 0:
        logging.warning("No %s ions – snapshot skipped.", ion_label)
        return

    # phosphate centres for plotting
    phos_atoms = {lab: u.select_atoms(f"resname I* and name {lab}").atoms[0]
                  for lab in PHOS_LABELS}

    # reference atoms for the plane (any three IP6 atom names)
    reference_atoms = {}
    for name in reference_triplet:
        sel = u.select_atoms(f"resname I* and name {name}")
        if len(sel) != 1:
            raise ValueError(f"Expected exactly one IP6 atom named '{name}', found {len(sel)}")
        reference_atoms[name] = sel.atoms[0]
    A, B, C = (reference_atoms[n] for n in reference_triplet)

    # 3. build basis
    rA = A.position
    v1 = B.position - rA
    v2 = C.position - rA
    normal = np.cross(v1, v2)
    if np.linalg.norm(normal) < 1e-6:
        raise ValueError("Reference atoms are collinear – cannot define a unique plane")
    e1 = v1 / np.linalg.norm(v1)
    e2 = np.cross(normal, e1)
    e2 /= np.linalg.norm(e2)

    # 4. project coords
    def project(X: np.ndarray) -> np.ndarray:
        R = X - rA
        return np.vstack((R @ e1, R @ e2)).T

    na_proj   = project(na_atoms.positions)
    phos_proj = project(np.array([phos_atoms[label].position for label in PHOS_LABELS]))

    # compute ion–peripheral-O minimum distances for dashed links
    periph_ag = {
        lab: u.select_atoms("resname I* and name " + " ".join(PHOS_OXYGENS[lab]))
        for lab in PHOS_LABELS
    }
    N = len(na_atoms)
    min_periph = np.full((N, len(PHOS_LABELS)), np.inf, dtype=np.float32)
    for j, lab in enumerate(PHOS_LABELS):
        d = distance_array(na_atoms.positions,
                           periph_ag[lab].positions,
                           box=u.dimensions)
        min_periph[:, j] = d.min(axis=1)
    closest_idx  = min_periph.argmin(axis=1)
    closest_dist = min_periph[np.arange(N), closest_idx]
    r_cut_A = None if r_cut_nm is None else r_cut_nm * 10.0

    # 5. plot
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")

    cmap = plt.cm.get_cmap("tab10", len(PHOS_LABELS))
    # legend for phosphates
    handles = []
    for j, (x, y) in enumerate(phos_proj):
        sc = ax.scatter(x, y, s=90, marker="o", edgecolor="k",
                        facecolor=cmap(j), zorder=3, label=PHOS_LABELS_CORRECTED[j])
        handles.append(sc)
    ax.legend(handles=handles, loc="best", frameon=True, fontsize=10)

    # Na ions
    ax.scatter(na_proj[:, 0], na_proj[:, 1],
               s=20, marker=".", color="black", zorder=2)

    # dashed links
    for i, (x_i, y_i) in enumerate(na_proj):
        j = closest_idx[i]
        if r_cut_A is not None and closest_dist[i] > r_cut_A:
            continue
        x_ph, y_ph = phos_proj[j]
        ax.plot([x_i, x_ph], [y_i, y_ph],
                ls="--", lw=0.6, color="grey", zorder=1)

    # labels & title
    ax.set_xlabel("e₁ projection / Å")
    ax.set_ylabel("e₂ projection / Å")
    tri = "–".join(reference_triplet)
    ax.set_title(f"plane ({tri}) – frame {frame}")
    ax.grid(True, ls=":", lw=0.4)

    plt.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    logging.info("Plane-projected snapshot (%s) written to %s", tri, out_png)


# ┃   3-D NETWORK  – ORIGINAL COORDS, CAMERA = normal(P,P2,P4)    ┃
# ╚═══════════════════════════════════════════════════════════════╝
def plot_frame_network_3d_fixed_view(
        u: mda.Universe,
        frame: int,
        out_png: pathlib.Path,
        anchors: Tuple[str, str, str] = ("P", "P2", "P4"),
        r_cut_nm: float = RCUTOFF_NM,
        look_along: int = +1,          # +1 → camera along  +n,  -1 → along –n
        ion_name: str = "NA",
        ion_label: Optional[str] = None,
) -> None:
    """
    Draw a 3-D snapshot in laboratory coordinates but orient the Matplotlib
    camera so that it looks *perpendicularly* onto the plane defined by
    the three anchor phosphates (default: P, P2, P4).

    Parameters
    ----------
    u          : MDAnalysis.Universe
    frame      : int          Frame index (0 ≤ frame < n_frames).
    out_png    : pathlib.Path Output image path.
    anchors    : (str,str,str)  Three phosphate labels fixing the plane.
    r_cut_nm   : float | None   Show dashed link only if ion–P distance
                                < r_cut_nm (None → show all).
    look_along : {+1,-1}        Choose whether the camera points parallel
                                (+1) or antiparallel (-1) to the normal.
    """
    ion_label = ion_label or f"{ion_name.upper()}⁺"
    # ── 1. go to requested frame ─────────────────────────────────
    if not (0 <= frame < len(u.trajectory)):
        raise IndexError(f"Frame {frame} out of range.")
    u.trajectory[frame]

    # ── 2. selections ────────────────────────────────────────────
    na_atoms = u.select_atoms(f"name {ion_name}")
    if len(na_atoms) == 0:
        logging.warning("No %s ions — skipping 3-D plot.", ion_label)
        return

    phos_atoms = {lab: u.select_atoms(f"resname I* and name {lab}").atoms[0]
                  for lab in PHOS_LABELS}

    try:
        r_P  = phos_atoms[anchors[0]].position
        r_P2 = phos_atoms[anchors[1]].position
        r_P4 = phos_atoms[anchors[2]].position
    except KeyError as e:
        raise ValueError(f"Anchor label {e.args[0]} not found among phosphates.")

    # ── 3. normal vector of anchor plane ─────────────────────────
    v1 = r_P2 - r_P
    v2 = r_P4 - r_P
    n = np.cross(v1, v2)
    if np.linalg.norm(n) < 1e-6:
        raise ValueError("Anchor atoms are colinear – plane undefined.")
    n /= np.linalg.norm(n)
    n *= look_along              # flip if the user wants the opposite side

    # convert normal into spherical view angles for Matplotlib
    elev_deg = np.degrees(np.arcsin(n[2]))              # arcsin(z/|n|)
    azim_rad = np.arctan2(n[1], n[0])                  # atan2(y,x)
    azim_deg = np.degrees(azim_rad)

    # ── 4. coordinates (original) & PERIPHERAL-oxygen distances ────
    na_pos   = na_atoms.positions                       # (N,3) Å
    phos_pos = np.array([phos_atoms[label].position for label in PHOS_LABELS])

    # Build a dict AtomGroup → peripheral Oxygens (three each)
    periph_ag = {label: u.select_atoms(
                    "resname I* and name " + " ".join(PHOS_OXYGENS[label]))
                for label in PHOS_LABELS}

    # matrix of minimum Na–O distances  (N_Na, 6)
    min_periph = np.full((len(na_atoms), 6), np.inf, dtype=np.float32)
    for j, label in enumerate(PHOS_LABELS):
        d = distance_array(na_pos, periph_ag[label].positions, box=u.dimensions)
        min_periph[:, j] = d.min(axis=1)          # min over the three O atoms

    closest_idx  = min_periph.argmin(axis=1)                      # phosphate index
    closest_dist = min_periph[np.arange(len(na_atoms)), closest_idx]
    r_cut_A      = r_cut_nm * 10.0                                # nm → Å
    # -------------------------------------------------------------------------

    # ── 5. plotting ─────────────────────────────────────────────
    fig = plt.figure(figsize=(6, 6))
    ax: Axes3D = fig.add_subplot(projection="3d")

    # ---- make all three coordinate panes pure white -----------------
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((1.0, 1.0, 1.0, 1.0))   # RGBA
        axis._axinfo['tick']['inward_factor']  = 0.0
        axis._axinfo['tick']['outward_factor'] = 0.0
        axis._axinfo['grid']['linewidth'] = 0.3
    # ----------------------------------------------------------------

    cmap = plt.cm.get_cmap("tab10", 6)
    # phosphates
    legend_handles = []                 # collect one handle per phosphate
    for j, (x, y, z) in enumerate(phos_pos):
        sc = ax.scatter(x, y, z, s=90, marker="o",
                        edgecolor="k", facecolor=cmap(j), depthshade=False,
                        label=PHOS_LABELS_CORRECTED[j])      # ← label for legend
        legend_handles.append(sc)                  # keep the handle
    
    ax.legend(handles=legend_handles,
          loc="lower right", frameon=True, fontsize=12) # title="Phosphates"

    # Na⁺ ions
    ax.scatter(na_pos[:, 0], na_pos[:, 1], na_pos[:, 2],
               s=15, c="black", marker=".", depthshade=False)

    # dashed connectors — draw *only* if Na is coordinated via peripheral O
    for i, (x_i, y_i, z_i) in enumerate(na_pos):
        if closest_dist[i] > r_cut_A:
            continue                    # not coordinated → no line
        j        = closest_idx[i]       # phosphate to which it is bound
        x_p, y_p, z_p = phos_pos[j]
        ax.plot([x_i, x_p], [y_i, y_p], [z_i, z_p],
                ls="--", lw=0.6, c="grey")

    # axis limits & labels
    span = np.ptp(np.vstack((na_pos, phos_pos)), axis=0)
    centre = np.mean(np.vstack((na_pos, phos_pos)), axis=0)
    max_half = span.max() / 2
    for axis, c in zip("xyz", centre):
        getattr(ax, f"set_{axis}lim")(c - max_half, c + max_half)

    ax.set_xlabel("x / Å")
    ax.set_ylabel("y / Å")
    ax.set_zlabel("z / Å")

    # --- fixed camera orientation ---
    ax.view_init(elev=elev_deg, azim=azim_deg)



    #ax.set_title(f"3-D Na–phosphate network (frame {frame})")
    ax.set_title(f"frame: {frame}", fontsize=14, loc="left")
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    logging.info("Fixed-view 3-D snapshot written to %s  (elev %.1f°, azim %.1f°)",
                 out_png, elev_deg, azim_deg)
