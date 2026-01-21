#!/usr/bin/env python3
"""
Na_IP6_coordination.py (CLI)

Iterates micro-state folders under a project's process directory, builds a
Boolean coordination tensor (N_ion, 6, N_frames) indicating whether each cation
is closer than RCUTOFF to any peripheral oxygen of P1..P6, aggregates to
time-series counts, and produces plots.

CLI usage:
    python -m electrofit_analysis.cli.coordination.Na_IP6_coordination \
            --project /path/to/project [--subdir process] [--determine-global-y] \
            [--rdf-y-max 1800] [--plot-projection] [--ion-name K] \
            [--rdf-show-first-shell] [--rdf-cutoff-mode {fixed,first-shell}]

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
    plot_coordination_boxplot,
    coordination_boxplot_metrics,
    plot_frame_network_3d_fixed_view, 
    plot_frame_network_plane, 
    plot_rdf_periphO_Na,
    first_shell_end,
    )

import seaborn as sns
from electrofit_analysis.cli.common import (
    resolve_stage,
    normalize_micro_name,
    resolve_run_and_analyze_dirs,
)
sns.set_context("talk")

# ╔═══════════════════════════════════════════════════════════════╗
# ┃                         CONSTANTS                             ┃
# ╚═══════════════════════════════════════════════════════════════╝

PROCESS_DIR_NAME   = "process"
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
        r_cut_nm: float | None = None,
        rdf_y_max_global: float = None,
        plot_projection: bool = False,
        precomputed_rdf: dict[str, Tuple[np.ndarray, np.ndarray]] | None = None,
        produce_boxplot: bool = True,
        rdf_show_first_shell: bool = False,
        rdf_cutoff_mode: str = "fixed",
        coord_cutoff_source: str | None = None,
        ion_name: str = "NA",
        ion_label: str | None = None,
) -> None:
    """
    For a single trajectory/topology pair:
      1. Build Boolean coordination tensor  A[na, phosph, frame].
      2. Sum over Na dimension → counts_per_frame[phosph, frame].
      3. Plot time-series.
    """

    # ── Load MD trajectory ───────────────────────────────────
    ion_name = ion_name.upper()
    ion_label = ion_label or f"{ion_name}⁺"
    u = mda.Universe(top_file.as_posix(), traj_file.as_posix())
    na_atoms = u.select_atoms(f"name {ion_name}")
    n_na = len(na_atoms)
    if n_na == 0:
        logging.warning("No %s atoms found – skipped.", ion_label)
        return
    logging.info("Loaded %d %s atoms and %d frames", n_na, ion_label, len(u.trajectory))

    print("Unique atom names in IP6:")
    print(sorted(set(u.select_atoms("resname I*").names)))

    phos_ag = {p: u.select_atoms("resname I* and name " + " ".join(names))
               for p, names in PHOS_OXYGENS.items()}

    # Sanity-check selections
    for p, ag in phos_ag.items():
        if len(ag) != 3:
            raise ValueError(f"Expected 3 peripheral O atoms for {p}, got {len(ag)}")

    # ── RDF curves (used for plotting and RDF-derived cutoffs) ─────────────────
    r_max_nm = 1.2
    nbins = 240
    if precomputed_rdf is None:
        rdf_results: dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        for label, oxy_ag in phos_ag.items():
            rdf = InterRDF(oxy_ag, na_atoms, range=(0, r_max_nm * 10), nbins=nbins)
            rdf.run()
            rdf_results[label] = (rdf.bins / 10.0, rdf.rdf)  # Å → nm
        # InterRDF iterates over the trajectory; reset before the explicit loop.
        u.trajectory[0]
    else:
        rdf_results = precomputed_rdf

    shell_ends: list[float] = []
    for label in PHOS_LABELS:
        r_vals, g_raw = rdf_results[label]
        shell_end = first_shell_end(np.asarray(r_vals, dtype=float), np.asarray(g_raw, dtype=float) * 3.0)
        if np.isfinite(shell_end):
            shell_ends.append(float(shell_end))
    shell_mean = float(np.mean(shell_ends)) if shell_ends else float("nan")
    shell_std = float(np.std(shell_ends, ddof=1)) if len(shell_ends) > 1 else 0.0

    cutoff_source = coord_cutoff_source or "fixed (CLI)"
    if r_cut_nm is None:
        cutoff_source = "RDF mean first-shell end (per microstate)"
        if shell_ends:
            r_cut_nm = shell_mean
        else:
            cutoff_source = "fallback fixed (no RDF minimum detected)"
            r_cut_nm = RCUTOFF_NM
            logging.warning(
                "Could not detect first-shell end; falling back to cutoff %.3f nm",
                float(r_cut_nm),
            )

    logging.info("Coordination cutoff used for counts: %.3f nm (%s)", float(r_cut_nm), cutoff_source)
    if shell_ends:
        logging.info(
            "First-shell end summary (n=%d): mean=%.3f nm, std=%.3f nm",
            len(shell_ends),
            shell_mean,
            shell_std,
        )

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
        title    = f"{ion_label} coordination counts – {dest_dir.parent.name}",
        ion_label=ion_label,
    )

    # Log the boxplot statistics derived from the cutoff-based coordination counts.
    # (This is independent of RDF integration cutoffs.)
    try:
        stats = coordination_boxplot_metrics(counts_ts)["per_phosphate"]
        for label in PHOS_LABELS:
            s = stats[label]
            logging.info(
                "Coordination count stats %s: mean=%.3f, std=%.3f, median=%.3f, q1=%.3f, q3=%.3f, whiskers=[%.3f, %.3f], outliers=%d",
                label,
                s["mean"],
                s["std"],
                s["median"],
                s["q1"],
                s["q3"],
                s["whisker_low"],
                s["whisker_high"],
                s["n_outliers"],
            )
    except Exception:
        logging.exception("Failed to compute coordination boxplot statistics.")

    if produce_boxplot:
        boxplot_path = dest_dir / "NaP_coordination_boxplot.pdf"
        plot_coordination_boxplot(
            counts_ts=counts_ts,
            out_png=boxplot_path,
            showfliers=True,
            ion_label=ion_label,
            ylabel=f"{ion_label} count",
        )
    if plot_projection: 
        snapshot_png = dest_dir / "network_frame_3d.png"
        plot_frame_network_3d_fixed_view(
            u, frame=521, out_png=snapshot_png, r_cut_nm=r_cut_nm, ion_name=ion_name, ion_label=ion_label
        )
        snapshot_png = dest_dir / "network_frame.png"
        plot_frame_network_plane(
            u, frame=525, out_png=snapshot_png,
            reference_triplet=("C", "C2", "C4"), r_cut_nm=r_cut_nm,
            ion_name=ion_name, ion_label=ion_label,
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
            show_shell_cutoff=rdf_show_first_shell,
            cutoff_mode=rdf_cutoff_mode,
            r_cut_nm=r_cut_nm,
            rdf_results_override=rdf_results,
            ion_name=ion_name,
            ion_label=ion_label,
        )
    else:
        plot_rdf_periphO_Na(
            u,
            out_png=rdf_png,
            r_max=1.2,
            nbins=240,
            y_max_global=1800.0,
            show_shell_cutoff=rdf_show_first_shell,
            cutoff_mode=rdf_cutoff_mode,
            r_cut_nm=r_cut_nm,
            rdf_results_override=rdf_results,
            ion_name=ion_name,
            ion_label=ion_label,
        )

    logging.info("Finished %s", dest_dir.parent.name)

# ╔═══════════════════════════════════════════════════════════════╗
# ┃                 DIRECTORY TRAVERSAL (CLI entry)              ┃
# ╚═══════════════════════════════════════════════════════════════╝
def main(
    project_dir: str,
    subdir: str = PROCESS_DIR_NAME,
    stage: str = "final",
    only: set[str] | None = None,
    determine_global_y: bool = False,
    rdf_y_max: float | None = None,
    plot_projection: bool = False,
    rdf_data_path: str | None = None,
    rep: int | None = None,
    boxplot: bool = True,
    ion_name: str = "NA",
    rdf_show_first_shell: bool = False,
    rdf_cutoff_mode: str = "fixed",
    project_mode: str = "auto",
    coord_cutoff_nm: float | None = None,
    coord_cutoff_nm_global: bool = False,
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
    ion_name : str
        Atom name of the cation to analyse (e.g. "NA", "K").
    rdf_show_first_shell : bool
        If True, draw a dotted guide at the first-shell end (RDF minimum after the first peak).
    rdf_cutoff_mode : str
        Either "fixed" (integrate to r_cut_nm) or "first-shell" (integrate to the detected first-shell end).
    project_mode : str
        How to interpret project_dir: 'project' (single project), 'collection' (directory containing projects),
        or 'auto' (detect based on presence of a process directory).
    coord_cutoff_nm : float | None
        Fixed cutoff (nm) for the explicit coordination/count algorithm. When None (default),
        the cutoff is derived from the RDF first-shell end (per microstate by default).
    coord_cutoff_nm_global : bool
        If True, compute one RDF-derived cutoff across the scan scope (project or collection) and use it for
        all explicit coordination/count analyses.
    """
    run_dir_name, analyze_base = resolve_stage(stage)
    only_norm = {normalize_micro_name(x) for x in only} if only else None
    ion_name = ion_name.upper()
    ion_label = f"{ion_name}⁺"

    base_path = pathlib.Path(project_dir).resolve()

    def detect_project_roots() -> list[pathlib.Path]:
        mode = (project_mode or "auto").strip().lower()
        if mode == "project":
            return [base_path]
        if mode == "collection":
            return sorted([p for p in base_path.iterdir() if p.is_dir() and (p / subdir).is_dir()])
        # auto
        if (base_path / subdir).is_dir():
            return [base_path]
        children = [p for p in base_path.iterdir() if p.is_dir() and (p / subdir).is_dir()]
        if children:
            return sorted(children)
        raise FileNotFoundError(f"Could not find '{subdir}' under {base_path} (project) or any child (collection).")

    project_roots = detect_project_roots()
    logging.info("Coordination project roots (%s): %s", project_mode, ", ".join(str(p) for p in project_roots))
    if coord_cutoff_nm_global and coord_cutoff_nm is None and only_norm is None:
        logging.warning(
            "Global cutoff requested without -m/--molecule filter; this may mix different molecules under '%s'.",
            subdir,
        )

    def micro_key(project_root: pathlib.Path, microstate_dir: pathlib.Path) -> str:
        if project_root == base_path:
            return microstate_dir.name
        rel = project_root.relative_to(base_path).as_posix()
        return f"{rel}/{microstate_dir.name}"

    def iter_microstates():
        for proj in project_roots:
            process_dir = proj / subdir
            if not process_dir.is_dir():
                continue
            for microstate in sorted(process_dir.iterdir()):
                if not microstate.is_dir():
                    continue
                if only_norm and microstate.name not in only_norm:
                    continue

                run_dir, analyze_base_dir = resolve_run_and_analyze_dirs(
                    microstate, stage, run_dir_name, analyze_base, rep
                )
                if not run_dir.is_dir():
                    continue

                if stage.strip().lower() == "remd":
                    traj = run_dir / "remd_center.xtc"
                    top = run_dir / "remd.tpr"
                else:
                    traj = run_dir / "md_center.xtc"
                    top = run_dir / "md.tpr"

                if not traj.is_file() or not top.is_file():
                    continue

                dest_dir = analyze_base_dir / "NaP_coordination"
                yield proj, microstate, traj, top, dest_dir

    rdf_cache_arrays: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]] | None = None
    rdf_cache_path = pathlib.Path(rdf_data_path).resolve() if rdf_data_path else None

    if rdf_cache_path and not rdf_cache_path.exists() and not (determine_global_y and rdf_y_max is None):
        raise FileNotFoundError(f"RDF cache {rdf_cache_path} not found. Run with --determine-global-y to create it.")

    Y_GLOBAL_MAX: float | None = None
    global_cutoff_nm: float | None = None

    # ---------------------------------------------------------------------------
    # Attempt to reuse a cache if provided
    # ---------------------------------------------------------------------------
    if rdf_cache_path and rdf_cache_path.exists():
        with rdf_cache_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        meta = payload.get("meta", {})
        if meta.get("global_y") is not None:
            try:
                Y_GLOBAL_MAX = float(meta["global_y"])
            except Exception:
                Y_GLOBAL_MAX = None
        ion_cached = meta.get("ion_name")
        if ion_cached and ion_cached.lower() != ion_name.lower():
            logging.warning(
                "RDF cache %s was generated for ion '%s' but current ion is '%s'. Ignoring cache.",
                rdf_cache_path, ion_cached, ion_name
            )
            rdf_cache_arrays = None
            Y_GLOBAL_MAX = None
        else:
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
            cached_cutoff = meta.get("coord_cutoff_nm_global")
            if cached_cutoff is not None:
                try:
                    global_cutoff_nm = float(cached_cutoff)
                except Exception:
                    global_cutoff_nm = None
            msg = f"Loaded RDF cache from {rdf_cache_path} (global y-limit = {Y_GLOBAL_MAX if Y_GLOBAL_MAX is not None else float('nan'):.2f})."
            logging.info(msg)
            print(msg)
            # Keep determine_global_y as-is: if the user explicitly requests a scan, refresh it.

    need_scan = False
    if determine_global_y and rdf_y_max is None:
        need_scan = True
    if coord_cutoff_nm_global and coord_cutoff_nm is None and global_cutoff_nm is None:
        need_scan = True

    if need_scan:
        # -----------------------------------------------------------------------
        # 0.  parameters for the scan
        # -----------------------------------------------------------------------
        r_max_nm = 1.2
        nbins = 240
        if determine_global_y and rdf_y_max is None:
            Y_GLOBAL_MAX = 0.0
        rdf_cache_arrays = {}
        shell_ends_global: list[float] = []

        # -----------------------------------------------------------------------
        # 1.  loop over all micro-state directories
        # -----------------------------------------------------------------------
        for proj, microstate, traj, top, dest_dir in iter_microstates():
            dest_dir.mkdir(exist_ok=True)
            logfile = dest_dir / "NaP_coordination.log"
            setup_logging(logfile)

            u = mda.Universe(top, traj, topology_format="TPR")
            ion_ag = u.select_atoms(f"name {ion_name}")
            if len(ion_ag) == 0:
                logging.warning("No %s atoms in %s – skipped.", ion_label, microstate.name)
                continue

            micro_cache: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
            for label in PHOS_LABELS:
                oxy_ag = u.select_atoms("resname I* and name " + " ".join(PHOS_OXYGENS[label]))
                rdf = InterRDF(oxy_ag, ion_ag, range=(0, r_max_nm * 10), nbins=nbins)
                rdf.run()

                r_nm = rdf.bins / 10.0
                g_raw = rdf.rdf
                micro_cache[label] = (r_nm.copy(), g_raw.copy())

                if determine_global_y and rdf_y_max is None and Y_GLOBAL_MAX is not None:
                    peak_height = float((g_raw * 3).max())
                    Y_GLOBAL_MAX = max(Y_GLOBAL_MAX, peak_height)

                if coord_cutoff_nm_global and coord_cutoff_nm is None and global_cutoff_nm is None:
                    shell_end = first_shell_end(r_nm, g_raw * 3.0)
                    if np.isfinite(shell_end):
                        shell_ends_global.append(float(shell_end))

            rdf_cache_arrays[micro_key(proj, microstate)] = micro_cache

        if determine_global_y and rdf_y_max is None:
            if rdf_cache_arrays and Y_GLOBAL_MAX is not None:
                Y_GLOBAL_MAX *= 1.05
                msg = f"Global y-limit to use in all figures: {Y_GLOBAL_MAX:.2f}"
                logging.info(msg)
                print(msg)
            else:
                logging.warning("No microstates processed during global y scan; using default y-limit.")
                Y_GLOBAL_MAX = None

        if coord_cutoff_nm_global and coord_cutoff_nm is None and global_cutoff_nm is None:
            if shell_ends_global:
                global_cutoff_nm = float(np.mean(shell_ends_global))
                std_nm = float(np.std(shell_ends_global, ddof=1)) if len(shell_ends_global) > 1 else 0.0
                logging.info(
                    "Global coordination cutoff (RDF mean first-shell end, n=%d): %.3f ± %.3f nm",
                    len(shell_ends_global),
                    global_cutoff_nm,
                    std_nm,
                )
            else:
                logging.warning("Could not determine global RDF cutoff; will fall back to per-microstate cutoffs.")

        if rdf_cache_path:
            payload = {
                "meta": {
                    "version": 2,
                    "r_max": r_max_nm,
                    "nbins": nbins,
                    "global_y": Y_GLOBAL_MAX,
                    "ion_name": ion_name,
                    "project_mode": project_mode,
                    "subdir": subdir,
                    "stage": stage,
                    "rep": rep,
                    "coord_cutoff_nm_global": global_cutoff_nm,
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

    if Y_GLOBAL_MAX is None:
        if rdf_y_max is not None:
            Y_GLOBAL_MAX = rdf_y_max
            logging.info(f"Using user-provided RDF y-limit: {Y_GLOBAL_MAX:.2f}")
        else:
            Y_GLOBAL_MAX = 1800.0
            logging.warning("No global y-limit determined; falling back to default = 1800.0")

    for proj, microstate, traj, top, dest_dir in iter_microstates():
        dest_dir.mkdir(exist_ok=True)
        logfile = dest_dir / "NaP_coordination.log"
        setup_logging(logfile)

        effective_y_max = rdf_y_max if rdf_y_max is not None else Y_GLOBAL_MAX
        rdf_override = None
        if rdf_cache_arrays:
            key = micro_key(proj, microstate)
            if key in rdf_cache_arrays:
                rdf_override = rdf_cache_arrays[key]

        coord_cutoff_to_use: float | None
        coord_cutoff_source: str | None = None
        if coord_cutoff_nm is not None:
            coord_cutoff_to_use = float(coord_cutoff_nm)
            coord_cutoff_source = "fixed (CLI)"
        elif coord_cutoff_nm_global and global_cutoff_nm is not None:
            coord_cutoff_to_use = float(global_cutoff_nm)
            coord_cutoff_source = "RDF mean first-shell end (global)"
        else:
            coord_cutoff_to_use = None

        analyse_one_microstate(
            traj_file = traj,
            top_file  = top,
            dest_dir  = dest_dir,
            r_cut_nm  = coord_cutoff_to_use,
            rdf_y_max_global = effective_y_max,
            plot_projection = plot_projection,
            precomputed_rdf = rdf_override,
            produce_boxplot = boxplot,
            rdf_show_first_shell = rdf_show_first_shell or (rdf_cutoff_mode == "first-shell"),
            rdf_cutoff_mode = rdf_cutoff_mode,
            coord_cutoff_source = coord_cutoff_source,
            ion_name = ion_name,
            ion_label = ion_label,
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
    parser.add_argument(
        "--no-boxplot",
        dest="boxplot",
        action="store_false",
        help="Skip the per-microstate coordination boxplot.",
    )
    parser.add_argument(
        "--ion-name",
        default="NA",
        help="Atom name of the cation to analyze (e.g. NA, K).",
    )
    parser.add_argument(
        "--rdf-show-first-shell",
        action="store_true",
        help="Draw a dotted line at the RDF first-shell end (first minimum after first peak).",
    )
    parser.add_argument(
        "--rdf-cutoff-mode",
        choices=["fixed", "first-shell"],
        default="fixed",
        help="RDF cutoff/integration mode: fixed (uses r_cut_nm) or first-shell (uses RDF minimum).",
    )
    parser.add_argument(
        "--project-mode",
        choices=["auto", "project", "collection"],
        default="auto",
        help="Interpret --project as a single project or a collection of projects (default: auto).",
    )
    parser.add_argument(
        "--coord-cutoff-nm",
        type=float,
        default=None,
        help=(
            "Fixed cutoff in nm for the explicit coordination/count algorithm. "
            "If omitted, the cutoff is estimated from the RDF first-shell end."
        ),
    )
    parser.add_argument(
        "--coord-cutoff-nm-global",
        action="store_true",
        help="Use one global RDF-derived cutoff across the scan scope (project or collection).",
    )
    args = parser.parse_args()

    main(
        project_dir=os.path.abspath(args.project),
        subdir=args.subdir,
        determine_global_y=args.determine_global_y,
        rdf_y_max=args.rdf_y_max,
        plot_projection=args.plot_projection,
        rdf_data_path=args.rdf_data,
        boxplot=args.boxplot,
        ion_name=args.ion_name,
        rdf_show_first_shell=args.rdf_show_first_shell,
        rdf_cutoff_mode=args.rdf_cutoff_mode,
        project_mode=args.project_mode,
        coord_cutoff_nm=args.coord_cutoff_nm,
        coord_cutoff_nm_global=args.coord_cutoff_nm_global,
    )
