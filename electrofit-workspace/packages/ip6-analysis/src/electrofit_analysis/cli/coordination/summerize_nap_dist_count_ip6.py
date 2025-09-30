#!/usr/bin/env python3
"""
summarize_nap_dist_count.py

Walks through IP_* subdirectories under a derived root directory, reading
Na+–phosphate distance and ion count files from analyze_final_sim/NaP_dist_count/,
and (optionally) NaP_coordination/NaP_coordination_bool.npy, computing summary
statistics and producing summary plots grouped by binary patterns (5 ones, 4 ones, 3 ones).

Now supports selective generation via CLI flags:

  Metrics (mutually non-exclusive; if none given, all metrics are run):
    --distance        Run distance summaries (distances_NaP*.xvg)
    --count           Run count summaries (ion_count_P*.xvg)
    --coordination    Run coordination summaries (NaP_coordination_bool.npy)

  Plot types:
    --plots means,box,overlaid,stats          # for distance/count (default: all four)
    --coord-plots box,stats                   # for coordination (default: box only)

Outputs (depending on flags):
  - distance_5ones.pdf, distance_4ones.pdf, distance_3ones.pdf                (means; one subplot per pattern)
  - count_5ones.pdf,    count_4ones.pdf,    count_3ones.pdf                   (means; one subplot per pattern)
  - distance_*_box.pdf, count_*_box.pdf                                       (boxplots per pattern)
  - distance_*_overlaid.pdf, count_*_overlaid.pdf                             (means overlaid per group)
  - distance_stats.xvg, count_stats.xvg, coordcount_stats.xvg                 (tables)
  - coordcount_*_box.pdf                                                       (coordination boxplots)

CLI Usage:
    python summarize_nap_dist_count.py --project /path/to/project --out /path/to/output_dir
    python summarize_nap_dist_count.py -p PROJ -o OUT --coordination
    python summarize_nap_dist_count.py -p PROJ -o OUT --distance --count --plots means,stats

Requires:
    matplotlib, numpy
"""

import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import to_rgba

# Define groups of patterns
GROUPS = {
    "5ones": ["011111", "101111", "111011", "111101"],
    "4ones": ["010111", "101101", "111001"],
    "3ones": ["010101", "101010", "101100", "111000"],
}

# Labels for the x-axis ticks:
P_LABELS = [f"P{i}" for i in range(1, 7)]


def load_mean_value(xvg_path):
    """Load an XVG file, return the mean of the second column (or NaN)."""
    try:
        data = np.loadtxt(xvg_path, comments=("@", "#"))
        return np.nan if data.size == 0 else data[:, 1].mean()
    except Exception:
        return np.nan


def load_values(xvg_path):
    """Load an XVG file, return the full second column as 1D array (or empty)."""
    try:
        data = np.loadtxt(xvg_path, comments=("@", "#"))
        return data[:, 1] if data.size else np.array([])
    except Exception:
        return np.array([])


def load_bool_counts(npy_path):
    """
    Read NaP_coordination_bool.npy and return a list of six numpy arrays,
    each containing the per–frame Na+ count for one phosphate (P1…P6).
    File stores a Boolean tensor of shape (N_Na, 6, N_frames).
    """
    if not os.path.isfile(npy_path):
        return [np.array([]) for _ in range(6)]
    try:
        A = np.load(npy_path, mmap_mode="r")  # read-only mem-map
        counts_ts = A.sum(axis=0)  # (6, N_frames)
        return [counts_ts[i] for i in range(6)]
    except Exception:
        return [np.array([]) for _ in range(6)]


# --- NEW: write stats for coordination ---
def write_stats_coord(root_dir, output_dir, metric_name="coordcount"):
    """
    Compute mean and std for each pattern and phosphate from Boolean coordination
    tensor (NaP_coordination_bool.npy) and save to <output_dir>/<metric_name>_stats.xvg.
    """
    header = ["# pattern"] + [f"P{i}_{stat}" for i in range(1, 7) for stat in ("mean", "std")]
    lines = [" ".join(header)]

    for group, patterns in GROUPS.items():
        for pat in patterns:
            npy_path = os.path.join(
                root_dir, f"IP_{pat}", "analyze_final_sim", "NaP_coordination", "NaP_coordination_bool.npy"
            )
            arrs = load_bool_counts(npy_path)
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
    """Compute mean and std for each pattern and phosphate to a .xvg file."""
    lines = []
    header = ["# pattern"] + [f"P{i}_{stat}" for i in range(1, 7) for stat in ("mean", "std")]
    lines.append(" ".join(header))
    for group, patterns in GROUPS.items():
        for pat in patterns:
            stats = []
            for suf in suffixes:
                fname = f"{file_prefix}{suf}.xvg"
                fpath = os.path.join(root_dir, f"IP_{pat}", "analyze_final_sim", "NaP_dist_count", fname)
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
    """Means per phosphate; subplots share y-limits."""
    for group, patterns in GROUPS.items():
        # Load to compute global y-range
        all_means = []
        for pat in patterns:
            means = []
            for suf in suffixes:
                fpath = os.path.join(root_dir, f"IP_{pat}", "analyze_final_sim", "NaP_dist_count", f"{file_prefix}{suf}.xvg")
                means.append(load_mean_value(fpath))
            all_means.append(means)

        flat = [m for sub in all_means for m in sub if not np.isnan(m)]
        if not flat:
            print(f"[WARN] No data for {metric_name} in group {group}; skipping.")
            continue
        ymin, ymax = min(flat), max(flat)
        pad = 0.1 * (ymax - ymin if ymax > ymin else 1.0)
        ymin -= pad; ymax += pad

        n = len(patterns)
        fig, axes = plt.subplots(nrows=n, ncols=1, figsize=(6, 1.5 * n), sharex=True, sharey=True)
        if n == 1:
            axes = [axes]

        for ax, pat, means in zip(axes, patterns, all_means):
            ax.plot(range(1, 7), means, marker="o", linestyle="-.", c="k")
            ax.set_title(pat, fontsize=10)
            ax.grid(True, axis="y", linestyle="--", alpha=0.6)
            ax.set_ylim(ymin, ymax)

        axes[-1].set_xticks(range(1, 7))
        axes[-1].set_xticklabels(P_LABELS, fontsize=10)
        for ax in axes: ax.set_ylabel("")
        fig.supxlabel("Phosphate site", fontsize=12)
        fig.supylabel(ylabel, fontsize=14)

        plt.tight_layout()
        out_pdf = os.path.join(output_dir, f"{metric_name}_{group}.pdf")
        fig.savefig(out_pdf)
        plt.close(fig)
        print(f"Saved {out_pdf}")


def summarize_metric_overlaid(root_dir, output_dir, metric_name, file_prefix, suffixes, ylabel):
    """For each GROUP, one figure with P1–P6 on x and one marker per pattern."""
    import itertools
    import matplotlib.cm as cm

    for group, patterns in GROUPS.items():
        all_means = {}
        flat = []
        for pat in patterns:
            means = []
            for suf in suffixes:
                fpath = os.path.join(
                    root_dir, f"IP_{pat}", "analyze_final_sim", "NaP_dist_count", f"{file_prefix}{suf}.xvg"
                )
                m = load_mean_value(fpath)
                means.append(m)
            all_means[pat] = means
            flat.extend([m for m in means if not np.isnan(m)])

        if not flat:
            print(f"[WARN] No data for {metric_name} (overlaid) in group {group}; skipping.")
            continue

        ymin, ymax = min(flat), max(flat)
        pad = 0.1 * (ymax - ymin if ymax > ymin else 1.0)
        ymin -= pad; ymax += pad

        fig, ax = plt.subplots(figsize=(6, 4))
        cmap = cm.get_cmap("tab10")
        colors = itertools.cycle(cmap.colors)

        for pat, color in zip(patterns, colors):
            ax.plot(range(1, 7), all_means[pat], marker="D", linestyle="", label=pat, color=color)

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


def summarize_metric_boxplot(root_dir, output_dir, metric_name, file_prefix, suffixes, ylabel, showfliers=False):
    """Boxplots per pattern; one subplot per pattern; six boxes per plot."""
    for group, patterns in GROUPS.items():
        # Load all data
        all_data = []
        for pat in patterns:
            per_p = []
            for suf in suffixes:
                fp = os.path.join(root_dir, f"IP_{pat}", "analyze_final_sim", "NaP_dist_count", f"{file_prefix}{suf}.xvg")
                vals = load_values(fp)
                per_p.append(vals)
            all_data.append(per_p)

        # Whisker-based global y-limits
        whisk_lows, whisk_highs = [], []
        for per_p in all_data:
            for arr in per_p:
                if arr.size:
                    q1, q3 = np.percentile(arr, (25, 75))
                    iqr = q3 - q1
                    lw_cand = arr[arr >= (q1 - 1.5 * iqr)]
                    hw_cand = arr[arr <= (q3 + 1.5 * iqr)]
                    if lw_cand.size: whisk_lows.append(lw_cand.min())
                    if hw_cand.size: whisk_highs.append(hw_cand.max())

        if not whisk_lows or not whisk_highs:
            print(f"[WARN] No data for {metric_name} (boxplot) in group {group}; skipping.")
            continue

        ymin, ymax = float(np.min(whisk_lows)), float(np.max(whisk_highs))
        pad = 0.05 * (ymax - ymin if ymax > ymin else 1.0)
        ymin -= pad; ymax += pad

        n = len(patterns)
        fig, axes = plt.subplots(n, 1, figsize=(6, 1.5 * n), sharex=True, sharey=False)
        axes = [axes] if n == 1 else list(axes)

        for ax, pat, per_p in zip(axes, patterns, all_data):
            bp = ax.boxplot(
                per_p, positions=range(1, 7), widths=0.6,
                showfliers=showfliers, patch_artist=False,
                medianprops={"visible": False},
                showmeans=True,
                meanprops={"marker": "D", "markeredgecolor": "black", "markerfacecolor": "black", "markersize": 6},
            )
            ax.set_title(pat, fontsize=14)
            ax.set_ylim(ymin, ymax)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))

        axes[-1].set_xticks(range(1, 7))
        axes[-1].set_xticklabels(P_LABELS, fontsize=12)
        for ax in axes: ax.set_ylabel("")

        fig.supylabel(ylabel, fontsize=14)
        plt.tight_layout()
        out_pdf = os.path.join(output_dir, f"{metric_name}_{group}_box.pdf")
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved boxplot {out_pdf}")


def summarize_boxplot_generic(root_dir, output_dir, metric_name, ylabel, data_loader, showfliers=False):
    """Shared-y boxplot figures (one per GROUP) using an arbitrary loader."""
    os.makedirs(output_dir, exist_ok=True)

    for group, patterns in GROUPS.items():
        print(f"--- Processing group '{group}' with {len(patterns)} patterns ---")
        all_data = []
        for pat in patterns:
            arrs = data_loader(pat)
            if not isinstance(arrs, (list, tuple)) or len(arrs) != 6:
                raise ValueError(f"data_loader('{pat}') must return 6 arrays; got {type(arrs)}")
            all_data.append(arrs)

        whisk_lows, whisk_highs = [], []
        for per_p in all_data:
            for arr in per_p:
                if arr.size == 0: continue
                q1, q3 = np.percentile(arr, (25, 75))
                iqr = q3 - q1
                arr_f = arr.astype(float, copy=False)
                lw_cand = arr_f[arr_f >= q1 - 1.5 * iqr]
                hw_cand = arr_f[arr_f <= q3 + 1.5 * iqr]
                if lw_cand.size: whisk_lows.append(lw_cand.min())
                if hw_cand.size: whisk_highs.append(hw_cand.max())

        if not whisk_lows or not whisk_highs:
            print(f"  ▶ No data for group '{group}', skipping.")
            continue

        ymin, ymax = min(whisk_lows), max(whisk_highs)
        pad = 0.05 * (ymax - ymin if ymax > ymin else 1.0)
        ymin -= pad; ymax += pad
        
        # Set fixed y-limits for coordination counts (this is a global metric for the project: /home/johannal96/Publications.nobackup/2025/electrofit-ip6-paper-2025/ip6-run)
        if metric_name == "coordcount":
            ymin = -0.3
            ymax = 6.3

        print(f"  ▶ y-range = [{ymin:.3g}, {ymax:.3g}]")

        n = len(patterns)
        fig, axes = plt.subplots(n, 1, figsize=(6, 1.5 * n), sharex=False, sharey=False, gridspec_kw={"hspace": 0.6})
        axes = np.atleast_1d(axes)

        for idx, (ax, pat, per_p) in enumerate(zip(axes, patterns, all_data)):
            colors = ["blue" if c == "1" else "black" for c in pat]
            bp = ax.boxplot(
                per_p, positions=range(1, 7), widths=0.6, showfliers=showfliers,
                patch_artist=True, medianprops={"visible": False},
                showmeans=True,
                meanprops={"marker": "D", "markeredgecolor": "black", "markerfacecolor": "black", "markersize": 6},
            )
            for ibox, box in enumerate(bp["boxes"]):
                box.set_edgecolor("black")
                box.set_facecolor(to_rgba("blue", alpha=0.3) if colors[ibox] == "blue" else "white")

            ax.set_title(pat, fontsize=14)
            ax.set_ylim(ymin - 0.5, ymax)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            if idx < n - 1:
                ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
            else:
                ax.set_xticks(range(1, 7))
                ax.set_xticklabels(P_LABELS, fontsize=12)
            ax.set_ylabel("")

        try:
            fig.supylabel(ylabel, fontsize=14)
        except AttributeError:
            fig.text(0.04, 0.5, ylabel, va="center", ha="center", rotation="vertical", fontsize=14)

        plt.tight_layout(rect=[0.05, 0, 1, 1])
        out_pdf = os.path.join(output_dir, f"{metric_name}_{group}_box.pdf")
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Saved boxplot → {out_pdf}")


# ------------------------ NEW: Driver helpers for flags ------------------------

def _parse_list(opt_str, allowed, default):
    """Parse a comma-separated option list; validate against allowed set."""
    if opt_str is None:
        return set(default)
    items = {s.strip().lower() for s in opt_str.split(",") if s.strip()}
    bad = items - set(allowed)
    if bad:
        raise ValueError(f"Unknown items {sorted(bad)}; allowed: {sorted(allowed)}")
    return items


def run_distance(root_dir, output_dir, plots):
    """Run selected plot types for distance metric."""
    suf = [str(i) for i in range(1, 7)]
    if "means" in plots:
        summarize_metric(root_dir, output_dir, "distance", "distances_NaP", suf, "Mean Na⁺–P distance (nm)")
    if "box" in plots:
        summarize_metric_boxplot(root_dir, output_dir, "distance", "distances_NaP", suf, "Na⁺–P distance (nm)")
    if "overlaid" in plots:
        summarize_metric_overlaid(root_dir, output_dir, "distance", "distances_NaP", suf, "Mean Na⁺–P distance (nm)")
    if "stats" in plots:
        write_stats(root_dir, output_dir, "distance", "distances_NaP", suf)


def run_count(root_dir, output_dir, plots):
    """Run selected plot types for count metric."""
    suf = ["", "1", "2", "3", "4", "5"]  # ''→P1, '1'→P2, … '5'→P6
    if "means" in plots:
        summarize_metric(root_dir, output_dir, "count", "ion_count_P", suf, "Mean Na⁺ count")
    if "box" in plots:
        summarize_metric_boxplot(root_dir, output_dir, "count", "ion_count_P", suf, "Na⁺ count", showfliers=True)
    if "overlaid" in plots:
        summarize_metric_overlaid(root_dir, output_dir, "count", "ion_count_P", suf, "Mean Na⁺ count")
    if "stats" in plots:
        write_stats(root_dir, output_dir, "count", "ion_count_P", suf)


def run_coordination(root_dir, output_dir, plots):
    """Run selected plot types for coordination metric."""
    def _coord_loader(pattern, _root=root_dir):
        npy_path = os.path.join(_root, f"IP_{pattern}", "analyze_final_sim", "NaP_coordination", "NaP_coordination_bool.npy")
        if not os.path.isfile(npy_path):
            print(f"[WARN] Boolean tensor not found for pattern {pattern}: {npy_path}")
        arrs = load_bool_counts(npy_path)
        if all(a.size == 0 for a in arrs):
            print(f"[WARN] EMPTY coord data for {pattern}: {npy_path}")
        return arrs

    if "box" in plots:
        summarize_boxplot_generic(root_dir, output_dir, metric_name="coordcount", ylabel="Na⁺ count (coordination)", data_loader=_coord_loader, showfliers=True)
    if "stats" in plots:
        write_stats_coord(root_dir, output_dir, metric_name="coordcount")


# ------------------------ Main -------------------------------------------------

def main(root_dir, output_dir, do_distance, do_count, do_coord, plots, coord_plots):
    os.makedirs(output_dir, exist_ok=True)

    # If none selected, run all (backward-compatible default)
    if not any([do_distance, do_count, do_coord]):
        do_distance = do_count = do_coord = True

    if do_distance:
        run_distance(root_dir, output_dir, plots)
    if do_count:
        run_count(root_dir, output_dir, plots)
    if do_coord:
        run_coordination(root_dir, output_dir, coord_plots)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Summarize Na–P distances, counts, and coordination across IP_* patterns."
    )
    parser.add_argument(
        "-p", "--project", required=True,
        help="Path to the project directory; root_dir is computed as PROJECT/process."
    )
    parser.add_argument(
        "-o", "--out", "--output", dest="output_dir", required=True,
        help="Output directory for summary PDFs and stats tables."
    )

    # Metric selectors: if none given, all are run
    parser.add_argument("--distance", action="store_true", help="Generate distance summaries only (unless combined with others).")
    parser.add_argument("--count", action="store_true", help="Generate ion-count summaries only (unless combined with others).")
    parser.add_argument("--coordination", action="store_true", help="Generate coordination summaries only (unless combined with others). By default generates boxplots.")

    # Plot-type toggles
    parser.add_argument(
        "--plots", default="means,box,overlaid,stats",
        help="Comma-separated list of plot types for distance/count metrics. Choices: means,box,overlaid,stats. Default: all."
    )
    parser.add_argument(
        "--coord-plots", default="box",
        help="Comma-separated list of plot types for coordination. Choices: box,stats. Default: box."
    )

    args = parser.parse_args()

    project_dir = os.path.abspath(args.project)
    root_dir = os.path.join(project_dir, "process")
    output_dir = os.path.abspath(args.output_dir)

    allowed = {"means", "box", "overlaid", "stats"}
    coord_allowed = {"box", "stats"}
    try:
        plots = _parse_list(args.plots, allowed, default=allowed)
        coord_plots = _parse_list(args.coord_plots, coord_allowed, default={"box"})
    except ValueError as e:
        raise SystemExit(f"Argument error: {e}")

    main(
        root_dir=root_dir,
        output_dir=output_dir,
        do_distance=args.distance,
        do_count=args.count,
        do_coord=args.coordination,
        plots=plots,
        coord_plots=coord_plots,
    )