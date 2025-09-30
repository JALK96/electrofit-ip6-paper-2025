#!/usr/bin/env python3
"""
Na_IP6_coordination.py (CLI)

Iterates micro-state folders under a project's process directory, builds a
Boolean coordination tensor (N_Na, 6, N_frames) indicating whether each Na⁺ is
closer than RCUTOFF to any peripheral oxygen of P1..P6, aggregates to
time-series counts, and produces plots.

CLI usage:
    python -m electrofit_analysis.cli.coordination.Na_IP6_coordination \
            --project /path/to/project [--subdir process] [--determine-global-y] \
            [--rdf-y-max 1800] [--plot-projection]

Outputs are written per microstate into analyze_final_sim/NaP_coordination/.

Date  : 2025-07-30
"""

from __future__ import annotations
import json
import os
import logging
import pathlib
from typing import Dict, Tuple

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array         #  (https://userguide.mdanalysis.org/examples/analysis/distances_and_contacts/distances_between_selections.html)
from MDAnalysis.analysis.rdf import InterRDF

from electrofit.infra.logging import setup_logging
from electrofit_analysis.viz.coord_helpers import (
    plot_counts_subplots, 
    plot_frame_network_3d_fixed_view, 
    plot_frame_network_plane, 
    plot_rdf_periphO_Na
    )

import seaborn as sns
sns.set_context("talk")

# ╔═══════════════════════════════════════════════════════════════╗
# ┃                         CONSTANTS                             ┃
# ╚═══════════════════════════════════════════════════════════════╝

PROCESS_DIR_NAME   = "process"
RUN_DIR_NAME       = "run_final_gmx_simulation"
PLOT_DIR_NAME      = "analyze_final_sim/NaP_coordination"
PLOT_NAME_TEMPLATE = "NaP_coordination_counts.png"

RCUTOFF_NM         = 0.32       # coordination threshold  (change if desired)
SAVE_3D_BOOL       = True       # turn off if disk space is an issue
BOOL_FILENAME      = "NaP_coordination_bool.npy"   # stored as numpy.memmap

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

# Set logging level
LOGLEVEL = logging.INFO

# ╔═══════════════════════════════════════════════════════════════╗
# ┃              CORE ANALYSIS ── ONE MICRO-STATE                ┃
# ╚═══════════════════════════════════════════════════════════════╝
def analyse_one_microstate(
        traj_file: pathlib.Path,
        top_file : pathlib.Path,
        dest_dir : pathlib.Path,
        r_cut_nm: float = RCUTOFF_NM,
        rdf_y_max_global: float = None,
        plot_projection: bool = False,
        precomputed_rdf: dict[str, Tuple[np.ndarray, np.ndarray]] | None = None,
) -> None:
    """
    For a single trajectory/topology pair:
      1. Build Boolean coordination tensor  A[na, phosph, frame].
      2. Sum over Na dimension → counts_per_frame[phosph, frame].
      3. Plot time-series.
    """

    # ── Load MD trajectory ───────────────────────────────────
    u = mda.Universe(top_file.as_posix(), traj_file.as_posix())
    na_atoms = u.select_atoms("name NA")
    n_na = len(na_atoms)
    if n_na == 0:
        logging.warning("No Na⁺ atoms found – skipped.")
        return
    logging.info("Loaded %d Na⁺ atoms and %d frames", n_na, len(u.trajectory))

    print("Unique atom names in IP6:")
    print(sorted(set(u.select_atoms("resname I*").names)))

    phos_ag = {p: u.select_atoms("resname I* and name " + " ".join(names))
               for p, names in PHOS_OXYGENS.items()}

    # Sanity-check selections
    for p, ag in phos_ag.items():
        if len(ag) != 3:
            raise ValueError(f"Expected 3 peripheral O atoms for {p}, got {len(ag)}")

    # ── Prepare Boolean memmap (optionally) ─────────────────
    n_frames   = len(u.trajectory)
    n_phos     = len(PHOS_LABELS)
    bool_shape = (n_na, n_phos, n_frames)          # (i, j, k)
    coord_mem  = None
    if SAVE_3D_BOOL:
        mmap_path = dest_dir / BOOL_FILENAME
        # open_memmap writes a valid .npy header; np.memmap does NOT
        coord_mem = np.lib.format.open_memmap(
            mmap_path, mode='w+', dtype=np.bool_, shape=bool_shape
        )
        coord_mem[:] = False      # initialise so file is fully allocated
        coord_mem.flush()
        logging.info(
            "Boolean tensor mapped to %s  (≈ %.1f MB).",
            mmap_path, coord_mem.nbytes / 1e6
        )

    # ── Frame-by-frame loop (vectorised in C under the hood) ──
    counts_ts = np.zeros((n_phos, n_frames), dtype=np.int16)  # for final plot

    for k, ts in enumerate(u.trajectory):
        na_pos = na_atoms.positions               # (N_Na, 3)
        min_dists = np.full((n_na, n_phos), np.inf, dtype=np.float32)

        # Compute per-phosphate minimum distances Naᵢ ↔ O_peripheral
        for j, p in enumerate(PHOS_LABELS):
            d = distance_array(na_pos, phos_ag[p].positions, box=ts.dimensions)  # (N_Na, 3)
            min_dists[:, j] = d.min(axis=1)   # (N_Na,)

        # For each Na, find closest phosphate (index of min along axis=1)
        closest_phos_idx = np.argmin(min_dists, axis=1)           # (N_Na,)
        closest_dist     = min_dists[np.arange(n_na), closest_phos_idx]

        # Boolean mask: coordinated if distance < cutoff
        coordin_mask = closest_dist < (r_cut_nm * 10)   # MDAnalysis uses Å
        # Note: *10 converts nm → Å, because Universe positions are in Å.

        # Update Boolean tensor and counts
        if SAVE_3D_BOOL:
            coord_mem[np.arange(n_na), closest_phos_idx, k] = coordin_mask

        # Aggregate counts per phosphate for this frame
        for j in range(n_phos):
            counts_ts[j, k] = np.count_nonzero(coordin_mask & (closest_phos_idx == j))

    # Flush memmap to disk
    if SAVE_3D_BOOL:
        coord_mem.flush()

    # ── Plot counts vs. time ─────────────────────────────────
    plot_counts_subplots(
        counts_ts=counts_ts,
        timestep_ps=u.trajectory.dt,
        out_png  = dest_dir / PLOT_NAME_TEMPLATE,
        title    = f"Na⁺ coordination counts – {dest_dir.parent.name}",
    )
    if plot_projection: 
        snapshot_png = dest_dir / "network_frame_3d.png"
        plot_frame_network_3d_fixed_view(u, frame=521, out_png=snapshot_png, r_cut_nm=r_cut_nm)
        snapshot_png = dest_dir / "network_frame.png"
        plot_frame_network_plane(
            u, frame=525, out_png=snapshot_png,
            reference_triplet=("C", "C2", "C4"), r_cut_nm=r_cut_nm
        )

    rdf_png = dest_dir / "rdf_Na_periphO.pdf"
    # global y max evaluated once (see outcommentd code below in main function)
    if rdf_y_max_global is not None:
        plot_rdf_periphO_Na(
            u,
            out_png=rdf_png,
            r_max=1.2,
            nbins=240,
            y_max_global=rdf_y_max_global,
            rdf_results_override=precomputed_rdf,
        )
    else:
        plot_rdf_periphO_Na(
            u,
            out_png=rdf_png,
            r_max=1.2,
            nbins=240,
            y_max_global=1800.0,
            rdf_results_override=precomputed_rdf,
        )

    logging.info("Finished %s", dest_dir.parent.name)

# ╔═══════════════════════════════════════════════════════════════╗
# ┃                 DIRECTORY TRAVERSAL (CLI entry)              ┃
# ╚═══════════════════════════════════════════════════════════════╝
def main(
    project_dir: str,
    subdir: str = PROCESS_DIR_NAME,
    determine_global_y: bool = False,
    rdf_y_max: float | None = None,
    plot_projection: bool = False,
    rdf_data_path: str | None = None,
) -> None:
    """Run coordination analysis across micro-states.

    Parameters
    ----------
    project_dir : str
        Path to the project root directory.
    subdir : str
        Subdirectory under project root to traverse (default: "process").
    determine_global_y : bool
        If True, scan all microstates to determine a global RDF y-limit.
    rdf_y_max : float | None
        If provided, use this as the RDF y-limit (overrides scanning).
    plot_projection : bool
        If True, also generate 2D/3D network projection snapshots.
    rdf_data_path : str | None
        Optional cache file for RDF data. When supplied, the command can reuse
        previously computed RDF curves (and the corresponding global y-limit)
        instead of recomputing them from scratch. When combined with
        ``--determine-global-y`` the cache is refreshed/created.
    """
    project_path = pathlib.Path(project_dir).resolve()
    process_dir  = project_path / subdir
    print(f"Process Dir: {process_dir}")

    rdf_cache_arrays: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]] | None = None
    rdf_cache_path = pathlib.Path(rdf_data_path).resolve() if rdf_data_path else None

    if rdf_cache_path and not rdf_cache_path.exists() and not (determine_global_y and rdf_y_max is None):
        raise FileNotFoundError(f"RDF cache {rdf_cache_path} not found. Run with --determine-global-y to create it.")

    Y_GLOBAL_MAX: float | None = None

    # ---------------------------------------------------------------------------
    # Attempt to reuse a cache if provided
    # ---------------------------------------------------------------------------
    if rdf_cache_path and rdf_cache_path.exists():
        with rdf_cache_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        meta = payload.get("meta", {})
        try:
            Y_GLOBAL_MAX = float(meta["global_y"])
        except KeyError as exc:
            raise KeyError(f"Missing 'global_y' in RDF cache {rdf_cache_path}") from exc
        cached_data = payload.get("data", {})
        rdf_cache_arrays = {
            micro: {
                label: (
                    np.asarray(entry["r"], dtype=float),
                    np.asarray(entry["g"], dtype=float),
                )
                for label, entry in micro_dict.items()
            }
            for micro, micro_dict in cached_data.items()
        }
        msg = f"Loaded RDF cache from {rdf_cache_path} (global y-limit = {Y_GLOBAL_MAX:.2f})."
        logging.info(msg)
        print(msg)
        determine_global_y = False  # cache already contains the scan results

    if determine_global_y and rdf_y_max is None:
        # -----------------------------------------------------------------------
        # 0.  parameters for the scan
        # -----------------------------------------------------------------------
        r_max_nm = 1.2
        nbins = 240
        Y_GLOBAL_MAX = 0.0
        rdf_cache_arrays = {}

        # -----------------------------------------------------------------------
        # 1.  loop over all micro-state directories
        # -----------------------------------------------------------------------
        for microstate in sorted(process_dir.iterdir()):
            if not microstate.is_dir():
                continue

            run_dir = microstate / RUN_DIR_NAME
            if not run_dir.is_dir():
                logging.debug("Skipping %s (no %s)", microstate.name, RUN_DIR_NAME)
                continue

            traj = run_dir / "md_center.xtc"
            top = run_dir / "md.tpr"
            if not traj.is_file() or not top.is_file():
                logging.warning("Missing files in %s – skipped.", run_dir)
                continue

            dest_dir = microstate / PLOT_DIR_NAME
            dest_dir.mkdir(exist_ok=True)
            logfile = dest_dir / "NaP_coordination.log"
            setup_logging(logfile)

            u = mda.Universe(top, traj, topology_format="TPR")
            na_ag = u.select_atoms("name NA")

            micro_cache: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
            for label in PHOS_LABELS:
                oxy_ag = u.select_atoms("resname I* and name " + " ".join(PHOS_OXYGENS[label]))
                rdf = InterRDF(na_ag, oxy_ag, range=(0, r_max_nm * 10), nbins=nbins)
                rdf.run()

                r_nm = rdf.bins / 10.0
                g_raw = rdf.rdf
                micro_cache[label] = (r_nm.copy(), g_raw.copy())

                peak_height = (g_raw * 3).max()
                Y_GLOBAL_MAX = max(Y_GLOBAL_MAX, float(peak_height))

            rdf_cache_arrays[microstate.name] = micro_cache

        if rdf_cache_arrays:
            Y_GLOBAL_MAX *= 1.05
            msg = f"Global y-limit to use in all figures: {Y_GLOBAL_MAX:.2f}"
            logging.info(msg)
            print(msg)

            if rdf_cache_path:
                payload = {
                    "meta": {
                        "version": 1,
                        "r_max": r_max_nm,
                        "nbins": nbins,
                        "global_y": Y_GLOBAL_MAX,
                    },
                    "data": {
                        micro: {
                            label: {
                                "r": r_vals.tolist(),
                                "g": g_vals.tolist(),
                            }
                            for label, (r_vals, g_vals) in micro_dict.items()
                        }
                        for micro, micro_dict in rdf_cache_arrays.items()
                    },
                }
                rdf_cache_path.parent.mkdir(parents=True, exist_ok=True)
                with rdf_cache_path.open("w", encoding="utf-8") as handle:
                    json.dump(payload, handle)
                logging.info("Saved RDF cache to %s", rdf_cache_path)
        else:
            logging.warning("No microstates processed during global y scan; using default y-limit.")
            Y_GLOBAL_MAX = None

    if Y_GLOBAL_MAX is None:
        if rdf_y_max is not None:
            Y_GLOBAL_MAX = rdf_y_max
            logging.info(f"Using user-provided RDF y-limit: {Y_GLOBAL_MAX:.2f}")
        else:
            Y_GLOBAL_MAX = 1800.0
            logging.warning("No global y-limit determined; falling back to default = 1800.0")

    for microstate in sorted(process_dir.iterdir()):
        if not microstate.is_dir():
            continue

        run_dir = microstate / RUN_DIR_NAME
        if not run_dir.is_dir():
            logging.debug("Skipping %s (no %s)", microstate.name, RUN_DIR_NAME)
            continue

        dest_dir = microstate / PLOT_DIR_NAME
        dest_dir.mkdir(exist_ok=True)

        traj = run_dir / "md_center.xtc"
        top  = run_dir / "md.tpr"

        if not traj.is_file() or not top.is_file():
            logging.warning("Missing files in %s – skipped.", run_dir)
            continue

        logfile = dest_dir / "NaP_coordination.log"
        setup_logging(logfile)

        effective_y_max = rdf_y_max if rdf_y_max is not None else Y_GLOBAL_MAX
        rdf_override = None
        if rdf_cache_arrays and microstate.name in rdf_cache_arrays:
            rdf_override = rdf_cache_arrays[microstate.name]
        analyse_one_microstate(
            traj_file = traj,
            top_file  = top,
            dest_dir  = dest_dir,
            r_cut_nm  = RCUTOFF_NM,
            rdf_y_max_global = effective_y_max,
            plot_projection = plot_projection,
            precomputed_rdf = rdf_override,
        )

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(
        description=(
            "Analyze Na–IP6 coordination across microstates under a project directory."
        )
    )
    parser.add_argument("-p", "--project", required=True, help="Path to the project root directory.")
    parser.add_argument(
        "--subdir",
        default=PROCESS_DIR_NAME,
        help='Subdirectory under project root to traverse (default: "process").',
    )
    parser.add_argument(
        "--determine-global-y",
        action="store_true",
        help="Scan all microstates to determine a common RDF y-limit.",
    )
    parser.add_argument(
        "--rdf-y-max",
        type=float,
        default=None,
        help="Override RDF y-limit with a fixed value (skips scanning).",
    )
    parser.add_argument(
        "--rdf-data",
        default=None,
        help="Path to an RDF cache file (JSON). If present, reuse cached curves;"
             " with --determine-global-y the cache is refreshed/created.",
    )
    parser.add_argument(
        "--plot-projection",
        action="store_true",
        help="Also generate 2D/3D network projection snapshots.",
    )
    args = parser.parse_args()

    main(
        os.path.abspath(args.project),
        args.subdir,
        args.determine_global_y,
        args.rdf_y_max,
        args.plot_projection,
        args.rdf_data,
    )
