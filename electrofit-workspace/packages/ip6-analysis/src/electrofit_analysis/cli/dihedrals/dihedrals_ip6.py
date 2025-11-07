#!/usr/bin/env python3
import os
import types
import logging
from pathlib import Path
from typing import Iterable, List, Tuple, Dict

import numpy as np
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals

from electrofit.infra.logging import setup_logging

from electrofit_analysis.viz.dihedral_helpers import (
    plot_time_series_summary_for_c_groups,
    plot_group_time_series,
    split_carbon_oxygen,
    plot_hist_offset_carbon,
    plot_kde_offset_carbon,
    plot_oxygen_chunks_hist,
    plot_oxygen_chunks_summary_hist,
    plot_oxygen_chunks_summary_kde
)

from electrofit_analysis.structure.util.common_util import ensure_dir
from electrofit_analysis.cli.common import resolve_stage, normalize_micro_name

# -----------------------
# Global style / config
# -----------------------
PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
sns.set_context("talk")
sns.set_style({"font.family": "serif", "font.serif": "Times New Roman"})

cf = types.SimpleNamespace()
cf.color_map = "viridis"
cf.dark_gray = "dimgray"



# -----------------------
# Domain: Discovery & IO
# -----------------------
def iter_molecule_process_dirs(project_path: Path, run_dir_name: str, analyze_base: str, only: set[str] | None = None) -> Iterable[Tuple[str, Path, Path]]:
    """
    Yield (species_id, run_final_dir, analyze_dir) for each molecule under project/process/*.
    """
    process_dir = project_path / "process"
    only_norm = {normalize_micro_name(x) for x in only} if only else None
    for entry in sorted(process_dir.iterdir()):
        if not entry.is_dir():
            continue
        if only_norm and entry.name not in only_norm:
            continue
        species_id = entry.name.replace("IP_", "")
        run_dir = entry / run_dir_name
        if not run_dir.is_dir():
            continue
        analyze_dir = entry / analyze_base / "dihedral"
        yield species_id, run_dir, analyze_dir


def load_universe(run_dir: Path) -> mda.Universe:
    gro = run_dir / "md.gro"
    xtc = run_dir / "md_center.xtc"
    if not gro.is_file() or not xtc.is_file():
        raise FileNotFoundError(f"Missing md.gro or md_center.xtc in {run_dir}")
    return mda.Universe(str(gro), str(xtc))


# -----------------------
# Domain: Dihedral setup
# -----------------------
def dihedral_groups_spec() -> List[Dict]:
    """Return the list of dihedral groups as in the original script."""
    return [
        {
            "name": "Ring Dihedrals",
            "dihedrals": [
                ("C", "C1", "C2", "C3"),
                ("C1", "C2", "C3", "C4"),
                ("C2", "C3", "C4", "C5"),
                ("C3", "C4", "C5", "C"),
                ("C4", "C5", "C", "C1"),
                ("C5", "C", "C1", "C2"),
            ],
        },
        {
            "name": "Dihedrals for C1",
            "dihedrals": [("C", "O", "P", "O6"), ("C", "O", "P", "O7"), ("C", "O", "P", "O8")],
        },
        {
            "name": "Dihedrals for C2",
            "dihedrals": [("C1", "O1", "P1", "O9"), ("C1", "O1", "P1", "O10"), ("C1", "O1", "P1", "O11")],
        },
        {
            "name": "Dihedrals for C3",
            "dihedrals": [("C2", "O2", "P2", "O12"), ("C2", "O2", "P2", "O13"), ("C2", "O2", "P2", "O14")],
        },
        {
            "name": "Dihedrals for C4",
            "dihedrals": [("C3", "O3", "P3", "O15"), ("C3", "O3", "P3", "O16"), ("C3", "O3", "P3", "O17")],
        },
        {
            "name": "Dihedrals for C5",
            "dihedrals": [("C4", "O4", "P4", "O18"), ("C4", "O4", "P4", "O19"), ("C4", "O4", "P4", "O20")],
        },
        {
            "name": "Dihedrals for C6",
            "dihedrals": [("C5", "O5", "P5", "O21"), ("C5", "O5", "P5", "O22"), ("C5", "O5", "P5", "O23")],
        },
    ]


def flatten_groups(groups: List[Dict]) -> Tuple[List[Tuple[str, str, str, str]], List[int], List[str]]:
    """
    Return:
      - dihedral_definitions: List of 4-tuples
      - dihedral_group_indices: index into 'groups' for each dihedral
      - group_names: names aligned to groups
    """
    dihedral_definitions: List[Tuple[str, str, str, str]] = []
    dihedral_group_indices: List[int] = []
    group_names: List[str] = []

    for g_idx, g in enumerate(groups):
        group_names.append(g["name"])
        for dih in g["dihedrals"]:
            dihedral_definitions.append(dih)
            dihedral_group_indices.append(g_idx)

    return dihedral_definitions, dihedral_group_indices, group_names


def select_unique_atom(u: mda.Universe, name: str) -> int:
    sel = u.select_atoms(f"name {name}")
    if len(sel) != 1:
        raise ValueError(f"Atom selection for '{name}' did not return exactly one atom (got {len(sel)}).")
    return sel[0].index


def build_dihedral_indices(u: mda.Universe, dihedral_definitions: List[Tuple[str, str, str, str]]) -> List[List[int]]:
    indices: List[List[int]] = []
    for dih in dihedral_definitions:
        ids = [select_unique_atom(u, atom_name) for atom_name in dih]
        indices.append(ids)
    return indices


# -----------------------
# Domain: Computation
# -----------------------
def compute_dihedral_angles(
    u: mda.Universe,
    dihedral_indices: List[List[int]],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns (times_ns, angles_deg) with angles shape (n_frames, n_dihedrals)
    """
    times_ps: List[float] = []
    angles_rad: List[List[float]] = []

    for ts in u.trajectory:
        pos = u.atoms.positions
        box = u.dimensions
        frame_angles: List[float] = []
        for A, B, C, D in dihedral_indices:
            phi = calc_dihedrals(
                pos[A][None, :],
                pos[B][None, :],
                pos[C][None, :],
                pos[D][None, :],
                box=box,
            )[0]
            frame_angles.append(float(phi))
        angles_rad.append(frame_angles)
        times_ps.append(float(ts.time))

    times_ns = np.asarray(times_ps) / 1000.0
    angles_deg = np.rad2deg(np.asarray(angles_rad))
    return times_ns, angles_deg


def group_dihedral_positions(dihedral_group_indices: List[int], n_groups: int) -> List[List[int]]:
    """Map group -> list of dihedral indices belonging to that group."""
    result: List[List[int]] = [[] for _ in range(n_groups)]
    for dih_idx, g_idx in enumerate(dihedral_group_indices):
        result[g_idx].append(dih_idx)
    return result


# -----------------------
# Orchestration
# -----------------------
def process_one_molecule(logger: logging.Logger, project_path: Path, species_id: str, run_dir: Path, outdir: Path):
    ensure_dir(outdir)
    u = load_universe(run_dir)

    groups = dihedral_groups_spec()
    dihedral_definitions, dihedral_group_indices, _group_names = flatten_groups(groups)
    dih_indices = build_dihedral_indices(u, dihedral_definitions)

    logger.info("Computing dihedral angles …")
    times_ns, angles_deg = compute_dihedral_angles(u, dih_indices)
    group_idx_map = group_dihedral_positions(dihedral_group_indices, len(groups))

    logger.info("Plotting time-series summaries …")
    plot_time_series_summary_for_c_groups(outdir, times_ns, angles_deg, groups, group_idx_map, dihedral_definitions)
    plot_group_time_series(outdir, times_ns, angles_deg, groups, group_idx_map, dihedral_definitions)

    logger.info("Preparing histogram/KDE plots …")
    carbon_data, oxygen_data = split_carbon_oxygen(groups, group_idx_map, dihedral_definitions, angles_deg)
    plot_hist_offset_carbon(outdir, carbon_data)
    plot_kde_offset_carbon(outdir, carbon_data, species_id)
    plot_oxygen_chunks_hist(outdir, oxygen_data, chunk_size=3)
    plot_oxygen_chunks_summary_hist(outdir, oxygen_data, chunk_size=3)
    plot_oxygen_chunks_summary_kde(outdir, oxygen_data, species_id, chunk_size=3)

    logger.info("Finished plots for %s", species_id)


def main(project_path_str: str | None = None, stage: str = 'final', only: set[str] | None = None):


    project_path = Path(project_path_str or PROJECT_PATH).resolve()

    run_dir_name, analyze_base = resolve_stage(stage)
    for species_id, run_dir, outdir in iter_molecule_process_dirs(project_path, run_dir_name, analyze_base, only=only):
        log_path = outdir / "dihedrals.log"
        setup_logging(log_path=log_path)
        logger = logging.getLogger(__name__)
        logger.info("Project: %s", project_path)
        logger.info("Processing %s", species_id)
        process_one_molecule(logger, project_path, species_id, run_dir, outdir)


if __name__ == "__main__":
    # Optional: parse --project from CLI, else use env/default
    import argparse

    p = argparse.ArgumentParser(description="Analyze dihedrals from final GMX runs.")
    p.add_argument("--project", type=str, default=None, help="Electrofit project root (contains 'process/').")
    args = p.parse_args()

    main(args.project)
