#!/usr/bin/env python3
"""
summarize_nap_dist_count.py

Walks through IP_* subdirectories under a derived root directory, reading
Na+–phosphate distance and ion count files from analyze_final_sim/NaP_dist_count/,
computing mean values per phosphate, and producing summary plots for three
groups of patterns (5 ones, 4 ones, 3 ones). Generates six vector PDFs:
  - distance_5ones.pdf, distance_4ones.pdf, distance_3ones.pdf
  - count_5ones.pdf,   count_4ones.pdf,   count_3ones.pdf
  - coordcount_5ones.pdf,  coordcount_4ones.pdf, coordcount_3ones.pdf

CLI Usage:
    python summerize_nap_dist_count.py --project /path/to/project --out /path/to/output_dir

Notes:
    The script expects the analysis inputs under ROOT = PROJECT/process, i.e.,
    it derives the "root" folder by appending "process" to the provided project
    directory. The output directory is provided explicitly and its structure is
    kept unchanged.

Requires:
    matplotlib, numpy
"""

import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import to_rgba

# sns.set_style("whitegrid", {
#    'font.family':'serif',
#    'font.serif':['Times New Roman']
# })

# sns.set_style("whitegrid")

# Define groups of patterns
GROUPS = {
    "5ones": ["011111", "101111", "111011", "111101"],
    "4ones": ["010111", "101101", "111001"],
    "3ones": ["010101", "101010", "101100", "111000"],
}

# Labels for the x-axis ticks:
P_LABELS = [f"P{i}" for i in range(1, 7)]


def load_mean_value(xvg_path):
    """
    Load an XVG file (whitespace delimited, comments start with @ or #),
    return the mean of the second column (or NaN if missing/empty).
    """
    try:
        data = np.loadtxt(xvg_path, comments=("@", "#"))
        return np.nan if data.size == 0 else data[:, 1].mean()
    except Exception:
        return np.nan


def load_values(xvg_path):
    """
    Load an XVG file and return the full array of its second column.
    Returns an empty array if file is missing or has no data.
    """
    try:
        data = np.loadtxt(xvg_path, comments=("@", "#"))
        return data[:, 1] if data.size else np.array([])
    except Exception:
        return np.array([])


def load_bool_counts(npy_path):
    """
    Read NaP_coordination_bool.npy and return a list of six numpy arrays,
    each containing the per–frame Na⁺‐count for one phosphate (P1…P6).
    The file stores a Boolean tensor of shape (N_Na, 6, N_frames).
    """
    if not os.path.isfile(npy_path):
        # keep interface: always return six arrays
        return [np.array([]) for _ in range(6)]
    try:
        A = np.load(npy_path, mmap_mode="r")  # read‑only mem‑map
        counts_ts = A.sum(axis=0)  # (6, N_frames)
        return [counts_ts[i] for i in range(6)]
    except Exception:
        return [np.array([]) for _ in range(6)]


# --- NEW HELPER for coordination stats ---
def write_stats_coord(root_dir, output_dir, metric_name="coordcount"):
    """
    Compute mean and std for each pattern and phosphate based on the
    Boolean coordination tensor (NaP_coordination_bool.npy) and save to
    <output_dir>/<metric_name>_stats.xvg

    The output header mirrors the format of write_stats() so downstream
    parsers remain compatible.
    """
    header = ["# pattern"] + [
        f"P{i}_{stat}" for i in range(1, 7) for stat in ("mean", "std")
    ]
    lines = [" ".join(header)]

    for group, patterns in GROUPS.items():
        for pat in patterns:
            npy_path_sub = os.path.join(
                root_dir,
                f"IP_{pat}",
                "analyze_final_sim",
                "NaP_coordination",
                "NaP_coordination_bool.npy",
            )
            npy_path_flat = os.path.join(
                root_dir, f"IP_{pat}", "analyze_final_sim", "NaP_coordination_bool.npy"
            )
            if os.path.isfile(npy_path_sub):
                arrs = load_bool_counts(npy_path_sub)
            else:
                arrs = load_bool_counts(npy_path_flat)

            stats = []
            for arr in arrs:
                if arr.size:
                    stats.append(f"{arr.mean():.5f}")
                    stats.append(f"{arr.std(ddof=1):.5f}")
                else:
                    stats.extend(["nan", "nan"])
            lines.append(" ".join([pat] + stats))

    out_file = os.path.join(output_dir, f"{metric_name}_stats.xvg")
    with open(out_file, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"Saved coordination stats file: {out_file}")


def write_stats(root_dir, output_dir, metric_name, file_prefix, suffixes):
    """
    Compute mean and std for each pattern and phosphate,
    and save to a text .xvg file.
    """
    lines = []
    # header: # pattern P1_mean P1_std ... P6_mean P6_std
    header = ["# pattern"] + [
        f"P{i}_{stat}" for i in range(1, 7) for stat in ("mean", "std")
    ]
    lines.append(" ".join(header))
    for group, patterns in GROUPS.items():
        for pat in patterns:
            stats = []
            for suf in suffixes:
                fname = f"{file_prefix}{suf}.xvg"
                fpath = os.path.join(
                    root_dir, f"IP_{pat}", "analyze_final_sim", "NaP_dist_count", fname
                )
                vals = load_values(fpath)
                if vals.size:
                    stats.append(f"{np.mean(vals):.5f}")
                    stats.append(f"{np.std(vals, ddof=1):.5f}")
                else:
                    stats.extend(["nan", "nan"])
            lines.append(" ".join([pat] + stats))
    out_file = os.path.join(output_dir, f"{metric_name}_stats.xvg")
    with open(out_file, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved stats file: {out_file}")


def summarize_metric(root_dir, output_dir, metric_name, file_prefix, suffixes, ylabel):
    """
    Same as before but now all subplots share y-limits.
    """
    for group, patterns in GROUPS.items():
        # First pass: load everything so we can find the global min/max
        all_means = []
        for pat in patterns:
            means = []
            for suf in suffixes:
                fname = f"{file_prefix}{suf}.xvg"
                fpath = os.path.join(
                    root_dir, f"IP_{pat}", "analyze_final_sim", "NaP_dist_count", fname
                )
                means.append(load_mean_value(fpath))
            all_means.append(means)

        # flatten and ignore NaNs
        flat = [m for sub in all_means for m in sub if not np.isnan(m)]
        ymin, ymax = min(flat), max(flat)

        # Optionally pad a little
        pad = 0.1 * (ymax - ymin)
        ymin -= pad
        ymax += pad

        # Create shared‐y subplots
        n = len(patterns)
        fig, axes = plt.subplots(
            nrows=n, ncols=1, figsize=(6, 1.5 * n), sharex=True, sharey=True
        )
        if n == 1:
            axes = [axes]

        # Plot each pattern
        for ax, pat, means in zip(axes, patterns, all_means):
            ax.plot(range(1, 7), means, marker="o", linestyle="-.", c="k")
            ax.set_title(pat, fontsize=10)
            ax.grid(True, axis="y", linestyle="--", alpha=0.6)
            ax.set_ylim(ymin, ymax)

        # only bottom axis gets the x‐ticks
        axes[-1].set_xticks(range(1, 7))
        axes[-1].set_xticklabels(P_LABELS, fontsize=10)

        # remove per‐axis y‐labels (we’ll add one global label)
        for ax in axes:
            ax.set_ylabel("")

        # global axis labels
        fig.supxlabel("Phosphate site", fontsize=12)
        fig.supylabel(ylabel, fontsize=14)

        plt.tight_layout()
        out_pdf = os.path.join(output_dir, f"{metric_name}_{group}.pdf")
        fig.savefig(out_pdf)
        plt.close(fig)
        print(f"Saved {out_pdf}")


def summarize_metric_overlaid(
    root_dir, output_dir, metric_name, file_prefix, suffixes, ylabel
):
    """
    For each GROUP, plot one figure with P1–P6 on x and one line per pattern,
    showing the mean value in different colors.
    """
    import itertools

    import matplotlib.cm as cm

    for group, patterns in GROUPS.items():
        # 1) load all means
        all_means = {}
        flat = []
        for pat in patterns:
            means = []
            for suf in suffixes:
                fpath = os.path.join(
                    root_dir,
                    f"IP_{pat}",
                    "analyze_final_sim",
                    "NaP_dist_count",
                    f"{file_prefix}{suf}.xvg",
                )
                m = load_mean_value(fpath)
                means.append(m)
            all_means[pat] = means
            flat.extend([m for m in means if not np.isnan(m)])

        # 2) compute common y‐limits
        ymin, ymax = min(flat), max(flat)
        pad = 0.1 * (ymax - ymin)
        ymin -= pad
        ymax += pad

        # 3) set up figure
        fig, ax = plt.subplots(figsize=(6, 4))
        cmap = cm.get_cmap("tab10")
        colors = itertools.cycle(cmap.colors)

        # 4) plot each pattern
        for pat, color in zip(patterns, colors):
            ax.plot(
                range(1, 7),
                all_means[pat],
                marker="D",
                linestyle="",
                label=pat,
                color=color,
            )

        ax.set_xticks(range(1, 7))
        ax.set_xticklabels(P_LABELS, fontsize=10)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel("Phosphate site", fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.6)
        ax.legend(title=group, fontsize=9, title_fontsize=10, loc="best")

        plt.tight_layout()
        out_pdf = os.path.join(output_dir, f"{metric_name}_{group}_overlaid.pdf")
        fig.savefig(out_pdf)
        plt.close(fig)
        print(f"Saved overlaid plot: {out_pdf}")


def summarize_metric_boxplot(
    root_dir, output_dir, metric_name, file_prefix, suffixes, ylabel, showfliers=False
):
    """
    Like summarize_metric, but draws boxplots instead of plotting means.
    Each subplot is a boxplot of the six phosphate distributions.
    """
    for group, patterns in GROUPS.items():
        # First pass: load all data into nested list [pattern][phosphate]
        all_data = []
        for pat in patterns:
            per_p = []
            for suf in suffixes:
                fn = f"{file_prefix}{suf}.xvg"
                fp = os.path.join(
                    root_dir, f"IP_{pat}", "analyze_final_sim", "NaP_dist_count", fn
                )
                vals = load_values(fp)
                per_p.append(vals)
            all_data.append(per_p)

        # Compute global min/max for y-limits (ignoring empty arrays)
        # Compute whisker endpoints per distribution
        whisker_lows = []
        whisker_highs = []
        for per_p in all_data:  # per pattern
            for arr in per_p:  # per phosphate
                if arr.size:
                    q1 = np.percentile(arr, 25)
                    q3 = np.percentile(arr, 75)
                    iqr = q3 - q1
                    # lower whisker: smallest point ≥ Q1 - 1.5*IQR
                    lw = arr[arr >= (q1 - 1.5 * iqr)].min()
                    # upper whisker: largest point ≤ Q3 + 1.5*IQR
                    hw = arr[arr <= (q3 + 1.5 * iqr)].max()
                    whisker_lows.append(lw)
                    whisker_highs.append(hw)

        # global whisker-based limits
        ymin = np.min(whisker_lows)
        ymax = np.max(whisker_highs)

        pad = 0.05 * (ymax - ymin)
        ymin -= pad
        ymax += pad

        # Create shared‐y boxplot subplots
        n = len(patterns)
        fig, axes = plt.subplots(n, 1, figsize=(6, 1.5 * n), sharex=True, sharey=False)
        if n == 1:
            axes = [axes]

        for ax, pat, per_p in zip(axes, patterns, all_data):
            # boxplot wants a sequence of arrays
            ax.boxplot(
                per_p,
                positions=range(1, 7),
                widths=0.6,
                showfliers=showfliers,
                patch_artist=False,
                medianprops={"visible": False},
                showmeans=True,
                meanprops={
                    "marker": "D",  # diamond marker
                    "markeredgecolor": "black",
                    "markerfacecolor": "black",
                    "markersize": 6,
                },
            )
            ax.set_title(pat, fontsize=14)
            ax.set_ylim(ymin, ymax)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3, prune=None))
            # ax.grid(axis='y', linestyle='--', alpha=0.6)
        for ax in fig.get_axes():
            # x‐axis tick labels
            ax.tick_params(axis="y", labelsize=12)
        # Bottom subplot: x-ticks
        axes[-1].set_xticks(range(1, 7))
        axes[-1].set_xticklabels(P_LABELS, fontsize=12)

        # Remove all per-axis y-labels
        for ax in axes:
            ax.set_ylabel("")

        # Add one global xlabel/ylabel
        # fig.supxlabel("Phosphate site", fontsize=12)
        fig.supylabel(ylabel, fontsize=14)

        plt.tight_layout()
        out_pdf = os.path.join(output_dir, f"{metric_name}_{group}_box.pdf")
        fig.savefig(out_pdf)
        plt.close(fig)
        print(f"Saved boxplot {out_pdf}")


def summarize_boxplot_generic(
    root_dir, output_dir, metric_name, ylabel, data_loader, showfliers=False
):
    """
    Build shared-y boxplot figures (one per GROUP) using arbitrary data.

    Parameters
    ----------
    root_dir : str
        Base directory (unused within function but kept for API compatibility).
    output_dir : str
        Where to save the PDF files.
    metric_name : str
        Prefix for output filenames.
    ylabel : str
        Global y-axis label.
    data_loader : callable
        data_loader(pattern: str) -> list of 6 numpy arrays
    showfliers : bool
        Whether to show outliers in the boxplot.
    """
    # 1) Ensure imports and globals are present
    try:
        GROUPS
        P_LABELS
    except NameError as e:
        raise NameError(
            "You must define GROUPS (dict of group→list of 6-digit patterns) "
            "and P_LABELS (list of 6 x-axis labels) in your module."
        ) from e

    # 2) Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    for group, patterns in GROUPS.items():
        print(f"--- Processing group '{group}' with {len(patterns)} patterns ---")
        # -------- load data --------
        all_data = []
        for pat in patterns:
            arrs = data_loader(pat)
            if not isinstance(arrs, (list, tuple)) or len(arrs) != 6:
                raise ValueError(
                    f"data_loader('{pat}') must return a list/tuple of 6 numpy arrays;"
                    f" got {type(arrs)} of length {len(arrs) if hasattr(arrs, '__len__') else 'unknown'}"
                )
            all_data.append(arrs)

        # -------- compute y-limits based on whiskers --------
        whisk_lows, whisk_highs = [], []
        for per_p in all_data:
            for arr in per_p:
                if arr is None or not hasattr(arr, "size"):
                    continue
                if arr.size == 0:
                    continue
                # compute 1st & 3rd quartiles
                q1, q3 = np.percentile(arr, (25, 75))
                iqr = q3 - q1
                arr_f = arr.astype(float, copy=False)
                # candidate points within the whisker range
                lw_cand = arr_f[arr_f >= q1 - 1.5 * iqr]
                hw_cand = arr_f[arr_f <= q3 + 1.5 * iqr]
                if lw_cand.size:
                    whisk_lows.append(lw_cand.min())
                if hw_cand.size:
                    whisk_highs.append(hw_cand.max())

        if not whisk_lows or not whisk_highs:
            print(f"  ▶ No data for group '{group}', skipping.")
            continue

        ymin, ymax = min(whisk_lows), max(whisk_highs)
        pad = 0.05 * (ymax - ymin)
        ymin -= pad
        ymax += pad
        print(f"  ▶ y-range = [{ymin:.3g}, {ymax:.3g}] (pad = {pad:.3g})")

        # -------- set up subplots --------
        n = len(patterns)
        fig, axes = plt.subplots(
            n,
            1,
            figsize=(6, 1.5 * n),
            sharex=False,
            sharey=False,
            gridspec_kw={"hspace": 0.6},  # ← increase this number to add more space
        )
        axes = np.atleast_1d(axes)

        # common spine styling
        for ax in axes:
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        for idx, (ax, pat, per_p) in enumerate(zip(axes, patterns, all_data)):
            # Validate pattern format
            if len(pat) != 6 or any(c not in "01" for c in pat):
                raise ValueError(
                    f"Pattern '{pat}' must be a 6-digit binary string (e.g. '010101')."
                )

            # determine edge colors
            colors = ["blue" if c == "1" else "black" for c in pat]

            # draw the boxplot
            bp = ax.boxplot(
                per_p,
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
            # color edges
            for ibox, box in enumerate(bp["boxes"]):
                # always opaque black edge
                box.set_edgecolor("black")

                if colors[ibox] == "blue":
                    # set the face to blue@30% opacity, edge stays at 100%
                    rgba_face = to_rgba("blue", alpha=0.3)
                    box.set_facecolor(rgba_face)
                else:
                    box.set_facecolor("white")

            ax.set_title(pat, fontsize=14)
            ax.set_ylim(ymin - 0.5, ymax)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.tick_params(axis="y", labelsize=12)

            # hide x-axis (ticks + spine) on all but last
            if idx < n - 1:
                ax.spines["bottom"].set_visible(True)
                ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
            else:
                # last subplot: show labels
                ax.spines["bottom"].set_visible(True)
                ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
                ax.set_xticks(range(1, 7))
                ax.set_xticklabels(P_LABELS, fontsize=12)

            # clear individual y-labels
            ax.set_ylabel("")

        # add a shared y-label
        try:
            fig.supylabel(ylabel, fontsize=14)
        except AttributeError:
            # for older matplotlib versions
            fig.text(
                0.04,
                0.5,
                ylabel,
                va="center",
                ha="center",
                rotation="vertical",
                fontsize=14,
            )

        plt.tight_layout(rect=[0.05, 0, 1, 1])  # leave space for ylabel

        # Save
        out_pdf = os.path.join(output_dir, f"{metric_name}_{group}_box.pdf")
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Saved boxplot → {out_pdf}")


def main(root_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # --- DISTANCES: files distances_NaP1.xvg … distances_NaP6.xvg
    summarize_metric(
        root_dir,
        output_dir,
        metric_name="distance",
        file_prefix="distances_NaP",
        suffixes=[str(i) for i in range(1, 7)],  # '1','2',…,'6'
        ylabel="Mean Na⁺–P distance (nm)",
    )

    # --- COUNTS: files ion_count_P.xvg, ion_count_P1.xvg … ion_count_P5.xvg
    summarize_metric(
        root_dir,
        output_dir,
        metric_name="count",
        file_prefix="ion_count_P",
        suffixes=["", "1", "2", "3", "4", "5"],  # ''→P1, '1'→P2, … '5'→P6
        ylabel="Mean Na⁺ count",
    )

    # Boxplots of distances
    summarize_metric_boxplot(
        root_dir,
        output_dir,
        metric_name="distance",
        file_prefix="distances_NaP",
        suffixes=[str(i) for i in range(1, 7)],
        ylabel="Na⁺–P distance (nm)",
    )
    # Boxplots of counts
    summarize_metric_boxplot(
        root_dir,
        output_dir,
        metric_name="count",
        file_prefix="ion_count_P",
        suffixes=["", "1", "2", "3", "4", "5"],
        ylabel="Na⁺ count",
        showfliers=True,
    )
    # Now write out the statistics tables:
    write_stats(
        root_dir,
        output_dir,
        metric_name="distance",
        file_prefix="distances_NaP",
        suffixes=[str(i) for i in range(1, 7)],
    )
    write_stats(
        root_dir,
        output_dir,
        metric_name="count",
        file_prefix="ion_count_P",
        suffixes=["", "1", "2", "3", "4", "5"],
    )

    # stats table for coordination counts
    write_stats_coord(root_dir, output_dir)

    # instead of summarize_metric(...)
    summarize_metric_overlaid(
        root_dir,
        output_dir,
        metric_name="distance",
        file_prefix="distances_NaP",
        suffixes=[str(i) for i in range(1, 7)],
        ylabel="Mean Na⁺–P distance (nm)",
    )
    # and likewise for counts
    summarize_metric_overlaid(
        root_dir,
        output_dir,
        metric_name="count",
        file_prefix="ion_count_P",
        suffixes=["", "1", "2", "3", "4", "5"],
        ylabel="Mean Na⁺ count",
    )

    # --- BOXPLOTS for Boolean coordination tensor -----------------
    def _coord_loader(pattern, _root=root_dir):
        """
        Return six arrays of Na⁺ counts for IP_<pattern>.

        Boolean tensor is expected at:
            <root_dir>/IP_<pattern>/analyze_final_sim/NaP_coordination/NaP_coordination_bool.npy
        """
        npy_path = os.path.join(
            _root,
            f"IP_{pattern}",
            "analyze_final_sim",
            "NaP_coordination",
            "NaP_coordination_bool.npy",
        )
        if not os.path.isfile(npy_path):
            print(f"[WARN] Boolean tensor not found for pattern {pattern}: {npy_path}")
        arrs = load_bool_counts(npy_path)
        if all(a.size == 0 for a in arrs):
            print(f"[WARN] EMPTY coord data for {pattern}: {npy_path}")
        return arrs

    print("Summarizing coordination counts as boxplots...")

    summarize_boxplot_generic(
        root_dir,
        output_dir,
        metric_name="coordcount",
        ylabel="Na⁺ count (coordination)",
        data_loader=_coord_loader,
        showfliers=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Summarize Na–P distances, counts, and coordination across IP_* patterns."
        )
    )
    parser.add_argument(
        "-p",
        "--project",
        required=True,
        help="Path to the project directory; root_dir is computed as PROJECT/process.",
    )
    parser.add_argument(
        "-o",
        "--out",
        "--output",
        dest="output_dir",
        required=True,
        help="Output directory for summary PDFs and stats tables.",
    )
    args = parser.parse_args()

    project_dir = os.path.abspath(args.project)
    root_dir = os.path.join(project_dir, "process")
    output_dir = os.path.abspath(args.output_dir)

    main(root_dir, output_dir)
