#!/usr/bin/env python3
from pathlib import Path
from typing import List, Tuple, Dict
import types

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde

from electrofit_analysis.structure.util.common_util import shift_atom_number, chunks


cf = types.SimpleNamespace()
cf.color_map = "viridis"
cf.dark_gray = "dimgray"

# -----------------------
# Plotting helpers
# -----------------------
def plot_time_series_summary_for_c_groups(
    outdir: Path,
    times_ns: np.ndarray,
    angles_deg: np.ndarray,
    dihedral_groups: List[Dict],
    group_dihedral_indices: List[List[int]],
    dihedral_definitions: List[Tuple[str, str, str, str]],
) -> None:
    groups_to_plot = [(i, g) for i, g in enumerate(dihedral_groups) if g["name"].startswith("Dihedrals for C")]
    if not groups_to_plot:
        return

    colors = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(groups_to_plot)))
    fig, axes = plt.subplots(nrows=len(groups_to_plot), ncols=1, figsize=(5, 2 * len(groups_to_plot)), sharex=True)
    if len(groups_to_plot) == 1:
        axes = [axes]

    for ax, (original_group_idx, group), color in zip(axes, groups_to_plot, colors):
        dihedral_idxs = group_dihedral_indices[original_group_idx]
        if not dihedral_idxs:
            continue
        first_idx = dihedral_idxs[0]
        dihedral_atoms = dihedral_definitions[first_idx]
        label = "-".join(shift_atom_number(a) for a in dihedral_atoms).replace(" ", "_")
        ax.plot(times_ns, angles_deg[:, first_idx], label=label, color=color, linestyle="", marker="o")
        ax.set_ylabel("Angle (°)", fontsize=18)
    axes[-1].set_xlabel("Time (ns)", fontsize=18)

    plt.tight_layout()
    plt.savefig(outdir / "dihedrals_summary.pdf")
    plt.close()


def plot_group_time_series(
    outdir: Path,
    times_ns: np.ndarray,
    angles_deg: np.ndarray,
    dihedral_groups: List[Dict],
    group_dihedral_indices: List[List[int]],
    dihedral_definitions: List[Tuple[str, str, str, str]],
) -> None:
    for g_idx, group in enumerate(dihedral_groups):
        dihedral_idxs = group_dihedral_indices[g_idx]
        if not dihedral_idxs:
            continue

        # C* groups: separate subplots for each dihedral
        if group["name"].startswith("Dihedrals for C"):
            n = len(dihedral_idxs)
            fig, axes = plt.subplots(nrows=n, ncols=1, figsize=(7, 2 * n), sharex=True)
            if n == 1:
                axes = [axes]
            colors_list = ["darkblue", "darkgreen", "darkred"]
            for i, dih_idx in enumerate(dihedral_idxs):
                ax = axes[i]
                label = "-".join(shift_atom_number(a) for a in dihedral_definitions[dih_idx]).replace(" ", "_")
                ax.plot(times_ns, angles_deg[:, dih_idx], label=label, color=colors_list[i % len(colors_list)],
                        linestyle="", marker="o")
                if i == 1:
                    ax.set_ylabel("Dihedral Angle (°)", fontsize=14)
                ax.legend(loc="best", fontsize=10)
            axes[-1].set_xlabel("Time (ns)", fontsize=14)
            plt.tight_layout()
            plt.savefig(outdir / f"dihedrals_{g_idx+1}_{group['name'].replace(' ', '_')}.pdf")
            plt.close()

        # Ring group: single figure with multiple lines
        else:
            plt.figure(figsize=(7, 4))
            colors_carbon = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, 6))
            for i, dih_idx in enumerate(dihedral_idxs):
                plt.plot(times_ns, angles_deg[:, dih_idx], label=f"τ {i+1}",
                         color=colors_carbon[i % len(colors_carbon)], alpha=1 - i / 10)
            plt.xlabel("Time (ns)", fontsize=14)
            plt.xlim((float(np.min(times_ns)), float(np.max(times_ns))))
            plt.ylabel("Dihedral Angle (°)", fontsize=14)
            plt.title(f"Dihedral Angles Over Time - {group['name']}", fontsize=16)
            plt.legend(loc="best", fontsize=10)
            plt.tight_layout()
            plt.savefig(outdir / f"dihedrals_{g_idx+1}_{group['name'].replace(' ', '_')}.pdf")
            plt.close()


def split_carbon_oxygen(
    dihedral_groups: List[Dict],
    group_dihedral_indices: List[List[int]],
    dihedral_definitions: List[Tuple[str, str, str, str]],
    angles_deg: np.ndarray,
) -> Tuple[List[Tuple[int, np.ndarray, str]], List[Tuple[int, np.ndarray, str]]]:
    carbon, oxygen = [], []
    for g_idx, group in enumerate(dihedral_groups):
        dihedral_idxs = group_dihedral_indices[g_idx]
        if group["name"].startswith("Ring"):
            target = carbon
        else:
            target = oxygen
        for d_idx in dihedral_idxs:
            label = "-".join(shift_atom_number(a) for a in dihedral_definitions[d_idx]).replace(" ", "_")
            target.append((d_idx, angles_deg[:, d_idx], label))
    return carbon, oxygen


def plot_hist_offset_carbon(outdir: Path, carbon_data: List[Tuple[int, np.ndarray, str]]) -> None:
    if not carbon_data:
        return
    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    colors = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(carbon_data)))
    base_lw, dec = 2.0, 0.25
    for i, (_, vals, label) in reversed(list(enumerate(carbon_data))):
        counts, edges = np.histogram(vals, bins=360, range=(-180, 180), density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        offset = 0.02 * i
        lw = max(base_lw - dec * i, 1.0)
        ax.plot(centers, counts + offset, drawstyle="steps-mid", label=label, color=colors[i], linewidth=lw)
    ax.set_yticks([])
    ax.set_ylabel("Density", fontsize=14)
    ax.set_xlabel("Angle (°)", fontsize=14)
    ax.set_xlim(-180, 180)
    ax.set_title("Carbon (Ring) Dihedral Hist with Vertical Offset", fontsize=16)
    ax.legend(fontsize=10)
    ax.grid(False)
    plt.savefig(outdir / "dihedral_hist_carbon_offset.pdf")
    plt.close()


def plot_kde_offset_carbon(outdir: Path, carbon_data: List[Tuple[int, np.ndarray, str]], species_id: str) -> None:
    if not carbon_data:
        return
    x_grid = np.linspace(-180, 180, 361)
    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    colors = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(carbon_data)))
    for i, (_, vals, label) in reversed(list(enumerate(carbon_data))):
        kde = gaussian_kde(vals)
        density = kde(x_grid)
        offset = 0.02 * i
        ax.plot(x_grid, density + offset, label=label, color=colors[i], linewidth=1.5)
    ax.set_yticks([])
    ax.set_ylabel("Density", fontsize=18)
    ax.set_xlabel("Angle (°)", fontsize=18)
    ax.set_xlim(-180, 180)
    ax.set_title(f"({species_id})", fontsize=18)
    ax.legend(fontsize=10)
    ax.grid(False)
    plt.savefig(outdir / "dihedral_kde_carbon_offset.pdf")
    plt.close()


def plot_oxygen_chunks_hist(outdir: Path, oxygen_data: List[Tuple[int, np.ndarray, str]], chunk_size: int = 3) -> None:
    if not oxygen_data:
        return
    for chunk_idx, chunk_data in enumerate(chunks(oxygen_data, chunk_size)):
        fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
        colors_list = ["darkblue", "darkgreen", "darkred"]
        base_lw, dec = 2.0, 0.5
        for i, (_, vals, label) in reversed(list(enumerate(chunk_data))):
            counts, edges = np.histogram(vals, bins=360, range=(-180, 180), density=True)
            centers = 0.5 * (edges[:-1] + edges[1:])
            offset = 0.01 * i
            lw = max(base_lw - dec * i, 1.0)
            ax.plot(centers, counts + offset, drawstyle="steps-mid", label=label, color=colors_list[i], linewidth=lw)
            ax.hlines(offset, xmin=-180, xmax=180, linestyles="--", colors="black", linewidth=1, alpha=1)
        ax.set_yticks([])
        ax.set_ylabel("Density", fontsize=14)
        ax.set_xlabel("Angle (°)", fontsize=14)
        ax.set_xlim(-180, 180)
        ax.set_title(f"Oxygen Dihedral Hist Offset - Group {chunk_idx + 1}", fontsize=16)
        ax.legend(fontsize=10)
        ax.grid(False)
        plt.savefig(outdir / f"dihedral_hist_oxygen_group{chunk_idx + 1}_offset.pdf")
        plt.close()


def plot_oxygen_chunks_summary_hist(outdir: Path, oxygen_data: List[Tuple[int, np.ndarray, str]], chunk_size: int = 3) -> None:
    if not oxygen_data:
        return
    oxygen_chunks = list(chunks(oxygen_data, chunk_size))
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    colors = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(oxygen_chunks)))
    for i, chunk in enumerate(oxygen_chunks):
        if not chunk:
            continue
        _, vals, label = chunk[0]
        counts, edges = np.histogram(vals, bins=360, range=(-180, 180), density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        offset = 0.01 * i
        ax.plot(centers, counts + offset, drawstyle="steps-mid", label=label, color=colors[i], linewidth=1.5)
    ax.set_yticks([])
    ax.set_xlabel("Angle (°)", fontsize=14)
    ax.set_ylabel("Density", fontsize=14)
    ax.set_xlim(-180, 180)
    ax.set_title("Summary: Oxygen Dihedrals", fontsize=16)
    ax.legend(fontsize=10)
    ax.grid(False)
    plt.savefig(outdir / "dihedral_hist_oxygen_summary.pdf")
    plt.close()


def plot_oxygen_chunks_summary_kde(
    outdir: Path, oxygen_data: List[Tuple[int, np.ndarray, str]], species_id: str, chunk_size: int = 3
) -> None:
    if not oxygen_data:
        return
    oxygen_chunks = list(chunks(oxygen_data, chunk_size))
    x_grid = np.linspace(-180, 180, 361)

    # KDE summary (no offset)
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    colors = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(oxygen_chunks)))
    for i, chunk in enumerate(oxygen_chunks):
        if not chunk:
            continue
        _, vals, label = chunk[0]
        sns.kdeplot(x=vals, bw_adjust=1, label=label, color=colors[i], linewidth=1.5, ax=ax, clip=(-180, 180))
    ax.set_yticks([])
    ax.set_xlabel("Angle (°)", fontsize=14)
    ax.set_ylabel("Density", fontsize=14)
    ax.set_xlim(-180, 180)
    ax.set_title("Summary: Oxygen Dihedrals (KDE)", fontsize=16)
    ax.legend(fontsize=10)
    ax.grid(False)
    plt.savefig(outdir / "dihedral_kde_oxygen_summary.pdf")
    plt.close()

    # KDE with vertical offsets
    fig2, ax2 = plt.subplots(figsize=(4, 4), constrained_layout=True)
    colors2 = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(oxygen_chunks)))
    for i, chunk in enumerate(oxygen_chunks):
        if not chunk:
            continue
        _, vals, label = chunk[0]
        kde = gaussian_kde(dataset=vals, bw_method="silverman")
        density = kde(x_grid)
        offset = 0.01 * i
        ax2.plot(x_grid, density + offset, label=label, color=colors2[i], linewidth=1.5)
    ax2.set_yticks([])
    ax2.set_xlabel("Angle (°)", fontsize=18)
    ax2.set_ylabel("Density", fontsize=18)
    ax2.set_xlim(-180, 180)
    ax2.set_title(f"({species_id})", fontsize=18)
    ax2.legend(fontsize=10)
    ax2.grid(False)
    plt.savefig(outdir / "dihedral_kde_offset_oxygen_summary.pdf")
    plt.close()