#!/usr/bin/env python3
"""
Na_IP6_coordination.py (CLI)

Iterates micro-state folders under a project's process directory and computes
coordination counts per phosphate using a *nearest-phosphate* assignment:

- For each ion and frame, compute the minimum distance to the 3 peripheral
  oxygens of each phosphate group.
- Assign the ion to the single phosphate with the smallest such distance,
  provided that this minimum distance is below the cutoff.

This produces a Boolean tensor (N_ion, 6, N_frames) that is aggregated into
time-series counts and plotted.

Optionally, an additional *non-exclusive contacts* tensor can be written, where
an ion can be in contact with multiple phosphates in the same frame (any
phosphate with min(peripheral-O) distance < cutoff).

CLI usage:
    python -m electrofit_analysis.cli.coordination.Na_IP6_coordination \
            --project /path/to/project [--subdir process] [--determine-global-y] \
            [--rdf-y-max 1800] [--plot-projection] [--ion-name K] \
            [--rdf-norm {rdf,density,none}] \
            [--rdf-show-first-shell] [--rdf-cutoff-mode {fixed,first-shell}] \
            [--rdf-scale {raw,per-phosphate}] \
            [--rdf-oxygen-mode {pooled,per-oxygen,both}]

Outputs are written per microstate into analyze_final_sim/IonP_coordination/.

Date  : 2025-07-30
"""

from __future__ import annotations
import json
import os
import logging
import pathlib
from typing import Dict, Tuple, Literal

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
PLOT_NAME_TEMPLATE = "IonP_coordination_counts.png"

RCUTOFF_NM         = 0.32       # coordination threshold  (change if desired)
SAVE_3D_BOOL       = True       # turn off if disk space is an issue
BOOL_FILENAME      = "IonP_coordination_bool.npy"   # legacy nearest-phosphate tensor (stored as numpy.memmap)
NEAREST_BOOL_ALIAS = "IonP_nearest_bool.npy"       # generic alias name for BOOL_FILENAME
CONTACTS_FILENAME  = "IonP_contacts_bool.npy"      # non-exclusive contacts tensor (stored as numpy.memmap)
RDF_CURVES_CSV     = "rdf_periphO_curves.csv"      # plotted RDF data points
RDF_SUMMARY_CSV    = "rdf_periphO_summary.csv"     # per-phosphate RDF/CN summary
RDF_CURVES_OXY_CSV = "rdf_periphO_curves_per_oxygen.csv"
RDF_SUMMARY_OXY_CSV = "rdf_periphO_summary_per_oxygen.csv"


def _rdf_scale_factor(rdf_scale: Literal["raw", "per-phosphate"], n_oxygen_sites: int = 3) -> float:
    if rdf_scale == "raw":
        return 1.0
    if rdf_scale == "per-phosphate":
        return float(n_oxygen_sites)
    raise ValueError(f"Unsupported rdf_scale '{rdf_scale}'")


def _rdf_scale_description(rdf_scale: Literal["raw", "per-phosphate"]) -> str:
    return "raw (per oxygen)" if rdf_scale == "raw" else "per-phosphate (raw × 3)"


def _rdf_plot_modes(rdf_oxygen_mode: Literal["pooled", "per-oxygen", "both"]) -> tuple[str, ...]:
    if rdf_oxygen_mode == "pooled":
        return ("pooled",)
    if rdf_oxygen_mode == "per-oxygen":
        return ("per-oxygen",)
    if rdf_oxygen_mode == "both":
        return ("pooled", "per-oxygen")
    raise ValueError(f"Unsupported rdf_oxygen_mode '{rdf_oxygen_mode}'")


def _link_or_copy_bool_tensor(src: pathlib.Path, dst: pathlib.Path) -> None:
    """Create dst as an alias of src (hardlink preferred, then symlink, then copy)."""
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        os.link(src, dst)
        return
    except Exception:
        pass
    try:
        # Relative symlink to keep the directory self-contained.
        dst.symlink_to(src.name)
        return
    except Exception:
        pass
    # Fallback: copy (file is typically small: N_ion*6*N_frames booleans)
    import shutil

    shutil.copyfile(src, dst)

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
        rdf_norm: Literal["rdf", "density", "none"] = "rdf",
        rdf_scale: Literal["raw", "per-phosphate"] = "per-phosphate",
        rdf_oxygen_mode: Literal["pooled", "per-oxygen", "both"] = "pooled",
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
    logging.info("RDF normalization mode: %s", rdf_norm)
    logging.info("RDF plotting scale: %s", _rdf_scale_description(rdf_scale))
    logging.info("RDF oxygen mode: %s", rdf_oxygen_mode)

    phos_ag = {p: u.select_atoms("resname I* and name " + " ".join(names))
               for p, names in PHOS_OXYGENS.items()}

    # Sanity-check selections
    for p, ag in phos_ag.items():
        if len(ag) != 3:
            raise ValueError(f"Expected 3 peripheral O atoms for {p}, got {len(ag)}")
        logging.info(
            "[RDF setup] %s: O names=%s | O indices=%s | O ids=%s",
            p,
            ",".join(PHOS_OXYGENS[p]),
            ",".join(map(str, ag.indices.tolist())),
            ",".join(map(str, ag.ids.tolist())),
        )

    # ── RDF curves (used for plotting and RDF-derived cutoffs) ─────────────────
    r_max_nm = 1.2
    nbins = 240
    mean_box_volume_nm3_from_rdf: float | None = None
    if precomputed_rdf is None:
        rdf_results: dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        for label, oxy_ag in phos_ag.items():
            rdf = InterRDF(oxy_ag, na_atoms, range=(0, r_max_nm * 10), nbins=nbins, norm=rdf_norm)
            rdf.run()
            rdf_results[label] = (rdf.results.bins / 10.0, rdf.results.rdf)  # Å → nm
            if mean_box_volume_nm3_from_rdf is None:
                vol_cum = getattr(rdf, "volume_cum", None)
                n_frames_rdf = getattr(rdf, "n_frames", None)
                if vol_cum is not None and n_frames_rdf:
                    mean_box_volume_nm3_from_rdf = float(vol_cum / n_frames_rdf) / 1000.0
            logging.info(
                "[RDF setup] Computing InterRDF for %s: n_O=%d vs n_%s=%d, range=[0, %.2f] nm, nbins=%d",
                label,
                len(oxy_ag),
                ion_name.upper(),
                len(na_atoms),
                r_max_nm,
                nbins,
            )
        # InterRDF iterates over the trajectory; reset before the explicit loop.
        u.trajectory[0]
    else:
        rdf_results = precomputed_rdf

    rdf_scale_factor = _rdf_scale_factor(rdf_scale)
    shell_ends: list[float] = []
    for label in PHOS_LABELS:
        r_vals, g_raw = rdf_results[label]
        shell_end = first_shell_end(
            np.asarray(r_vals, dtype=float),
            np.asarray(g_raw, dtype=float) * rdf_scale_factor,
        )
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
    contacts_mem = None
    if SAVE_3D_BOOL:
        mmap_path = dest_dir / BOOL_FILENAME
        contacts_path = dest_dir / CONTACTS_FILENAME
        # open_memmap writes a valid .npy header; np.memmap does NOT
        coord_mem = np.lib.format.open_memmap(
            mmap_path, mode='w+', dtype=np.bool_, shape=bool_shape
        )
        coord_mem[:] = False      # initialise so file is fully allocated
        coord_mem.flush()
        logging.info(
            "Nearest-phosphate tensor mapped to %s  (≈ %.1f MB).",
            mmap_path, coord_mem.nbytes / 1e6
        )
        contacts_mem = np.lib.format.open_memmap(
            contacts_path, mode="w+", dtype=np.bool_, shape=bool_shape
        )
        contacts_mem[:] = False
        contacts_mem.flush()
        logging.info(
            "Contacts tensor mapped to %s  (≈ %.1f MB).",
            contacts_path, contacts_mem.nbytes / 1e6
        )

    # ── Frame-by-frame loop (vectorised in C under the hood) ──
    counts_ts = np.zeros((n_phos, n_frames), dtype=np.int16)  # for final plot
    volume_sum_nm3 = 0.0
    volume_samples = 0

    for k, ts in enumerate(u.trajectory):
        na_pos = na_atoms.positions               # (N_Na, 3)
        min_dists = np.full((n_na, n_phos), np.inf, dtype=np.float32)
        volume_sum_nm3 += float(np.prod(ts.dimensions[:3]) / 1000.0)  # Å^3 -> nm^3
        volume_samples += 1

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
            # Non-exclusive contacts: any phosphate whose min(peripheral-O) distance is < cutoff
            contacts_mem[:, :, k] = min_dists < (r_cut_nm * 10)

        # Aggregate counts per phosphate for this frame
        for j in range(n_phos):
            counts_ts[j, k] = np.count_nonzero(coordin_mask & (closest_phos_idx == j))

    mean_box_volume_nm3 = (
        volume_sum_nm3 / volume_samples if volume_samples > 0 else mean_box_volume_nm3_from_rdf
    )
    if mean_box_volume_nm3 is None and mean_box_volume_nm3_from_rdf is not None:
        mean_box_volume_nm3 = mean_box_volume_nm3_from_rdf
    ion_density_nm3 = (
        float(n_na / mean_box_volume_nm3) if mean_box_volume_nm3 and mean_box_volume_nm3 > 0 else None
    )
    logging.info(
        "[RDF setup] mean_box_volume=%.3f nm^3, n_%s=%d, ion_density=%.6f nm^-3",
        float(mean_box_volume_nm3) if mean_box_volume_nm3 is not None else float("nan"),
        ion_name.upper(),
        n_na,
        float(ion_density_nm3) if ion_density_nm3 is not None else float("nan"),
    )

    # Flush memmap to disk
    if SAVE_3D_BOOL:
        coord_mem.flush()
        contacts_mem.flush()
        try:
            _link_or_copy_bool_tensor(dest_dir / BOOL_FILENAME, dest_dir / NEAREST_BOOL_ALIAS)
            logging.info("Wrote alias tensor name %s -> %s", NEAREST_BOOL_ALIAS, BOOL_FILENAME)
        except Exception:
            logging.exception("Failed to create alias %s for %s", NEAREST_BOOL_ALIAS, BOOL_FILENAME)

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
        boxplot_path = dest_dir / "IonP_coordination_boxplot.pdf"
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

    y_for_plot = rdf_y_max_global if rdf_y_max_global is not None else 1800.0
    for mode in _rdf_plot_modes(rdf_oxygen_mode):
        if mode == "pooled":
            rdf_png = dest_dir / f"rdf_{ion_name.upper()}_periphO.pdf"
            rdf_curve_csv = dest_dir / RDF_CURVES_CSV
            rdf_summary_csv = dest_dir / RDF_SUMMARY_CSV
        else:
            rdf_png = dest_dir / f"rdf_{ion_name.upper()}_periphO_per_oxygen.pdf"
            rdf_curve_csv = dest_dir / RDF_CURVES_OXY_CSV
            rdf_summary_csv = dest_dir / RDF_SUMMARY_OXY_CSV

        plot_rdf_periphO_Na(
            u,
            out_png=rdf_png,
            r_max=1.2,
            nbins=240,
            y_max_global=y_for_plot,
            show_shell_cutoff=rdf_show_first_shell,
            cutoff_mode=rdf_cutoff_mode,
            r_cut_nm=r_cut_nm,
            rdf_results_override=rdf_results,
            rdf_norm=rdf_norm,
            rdf_scale=rdf_scale,
            rdf_oxygen_mode=mode,
            ion_name=ion_name,
            ion_label=ion_label,
            ion_density_nm3=ion_density_nm3,
            mean_volume_nm3=mean_box_volume_nm3,
            rdf_curve_csv=rdf_curve_csv,
            rdf_summary_csv=rdf_summary_csv,
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
    rdf_norm: Literal["rdf", "density", "none"] = "rdf",
    rdf_scale: Literal["raw", "per-phosphate"] = "per-phosphate",
    rdf_oxygen_mode: Literal["pooled", "per-oxygen", "both"] = "pooled",
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
    rdf_norm : {"rdf", "density", "none"}
        InterRDF normalization mode:
        - rdf: normalized radial distribution function g(r)
        - density: single-particle density n(r)
        - none: raw shell counts
    rdf_scale : {"raw", "per-phosphate"}
        RDF scaling used for plotting/integration:
        - raw: use InterRDF output directly (per oxygen)
        - per-phosphate: multiply raw RDF by 3 (sum over three peripheral O sites)
    rdf_oxygen_mode : {"pooled", "per-oxygen", "both"}
        RDF grouping mode for visualization:
        - pooled: one curve per phosphate using all three peripheral oxygens
        - per-oxygen: three curves per phosphate (one curve per terminal oxygen)
        - both: produce both output variants
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
    rdf_norm = rdf_norm.lower()
    ion_label = f"{ion_name}⁺"
    scale_factor = _rdf_scale_factor(rdf_scale)
    logging.info("RDF normalization mode: %s", rdf_norm)
    logging.info("RDF plotting scale: %s", _rdf_scale_description(rdf_scale))
    logging.info("RDF oxygen mode: %s", rdf_oxygen_mode)

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

                dest_dir = analyze_base_dir / "IonP_coordination"
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
            cached_norm = meta.get("rdf_norm")
            if cached_norm and cached_norm != rdf_norm:
                logging.warning(
                    "RDF cache norm was '%s' but current run uses '%s'. Ignoring cache.",
                    cached_norm,
                    rdf_norm,
                )
                rdf_cache_arrays = None
                Y_GLOBAL_MAX = None
                cached_data = {}
            cached_scale = meta.get("rdf_scale")
            if cached_data and cached_scale and cached_scale != rdf_scale:
                logging.info(
                    "RDF cache scale was '%s' but current run uses '%s'; raw curves are reused and y-limit is re-derived.",
                    cached_scale,
                    rdf_scale,
                )
            if cached_data:
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
            if cached_data and cached_cutoff is not None:
                try:
                    global_cutoff_nm = float(cached_cutoff)
                except Exception:
                    global_cutoff_nm = None
            if cached_data:
                msg = f"Loaded RDF cache from {rdf_cache_path} (global y-limit = {Y_GLOBAL_MAX if Y_GLOBAL_MAX is not None else float('nan'):.2f})."
                logging.info(msg)
                print(msg)
            # Keep determine_global_y as-is: if the user explicitly requests a scan, refresh it.

    if rdf_cache_arrays and rdf_y_max is None and (not determine_global_y):
        peak_max = 0.0
        for micro_dict in rdf_cache_arrays.values():
            for label in PHOS_LABELS:
                _, g_raw = micro_dict[label]
                peak_max = max(peak_max, float(np.max(g_raw * scale_factor)))
        if peak_max > 0.0:
            Y_GLOBAL_MAX = peak_max * 1.05
            logging.info(
                "Global y-limit derived from cached RDFs with scale '%s': %.2f",
                rdf_scale,
                Y_GLOBAL_MAX,
            )

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
            logfile = dest_dir / "IonP_coordination.log"
            setup_logging(logfile)

            u = mda.Universe(top, traj, topology_format="TPR")
            ion_ag = u.select_atoms(f"name {ion_name}")
            if len(ion_ag) == 0:
                logging.warning("No %s atoms in %s – skipped.", ion_label, microstate.name)
                continue

            micro_cache: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
            for label in PHOS_LABELS:
                oxy_ag = u.select_atoms("resname I* and name " + " ".join(PHOS_OXYGENS[label]))
                rdf = InterRDF(oxy_ag, ion_ag, range=(0, r_max_nm * 10), nbins=nbins, norm=rdf_norm)
                rdf.run()

                r_nm = rdf.results.bins / 10.0
                g_raw = rdf.results.rdf
                micro_cache[label] = (r_nm.copy(), g_raw.copy())

                if determine_global_y and rdf_y_max is None and Y_GLOBAL_MAX is not None:
                    peak_height = float((g_raw * scale_factor).max())
                    if rdf_oxygen_mode in {"per-oxygen", "both"}:
                        for oxygen_index in oxy_ag.indices.tolist():
                            oxygen_ag = u.atoms[[int(oxygen_index)]]
                            rdf_single = InterRDF(oxygen_ag, ion_ag, range=(0, r_max_nm * 10), nbins=nbins, norm=rdf_norm)
                            rdf_single.run()
                            peak_height = max(peak_height, float((rdf_single.results.rdf * scale_factor).max()))
                    Y_GLOBAL_MAX = max(Y_GLOBAL_MAX, peak_height)

                if coord_cutoff_nm_global and coord_cutoff_nm is None and global_cutoff_nm is None:
                    shell_end = first_shell_end(r_nm, g_raw * scale_factor)
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
                    "rdf_norm": rdf_norm,
                    "rdf_scale": rdf_scale,
                    "rdf_oxygen_mode": rdf_oxygen_mode,
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
        logfile = dest_dir / "IonP_coordination.log"
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
            rdf_norm = rdf_norm,
            rdf_scale = rdf_scale,
            rdf_oxygen_mode = rdf_oxygen_mode,
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
        help="Atom name of the cation to analyze (e.g. NA, K, MG, CA).",
    )
    parser.add_argument(
        "--rdf-norm",
        choices=["rdf", "density", "none"],
        default="rdf",
        help="InterRDF normalization mode: rdf (g(r), default), density (n(r)), or none (raw shell counts).",
    )
    parser.add_argument(
        "--rdf-scale",
        choices=["raw", "per-phosphate"],
        default="per-phosphate",
        help=(
            "RDF scaling for plotting/integration: raw (per oxygen) or per-phosphate (raw×3)."
        ),
    )
    parser.add_argument(
        "--rdf-oxygen-mode",
        choices=["pooled", "per-oxygen", "both"],
        default="pooled",
        help=(
            "RDF oxygen grouping: pooled (default), per-oxygen (three curves per phosphate), or both."
        ),
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
        rdf_norm=args.rdf_norm,
        rdf_scale=args.rdf_scale,
        rdf_oxygen_mode=args.rdf_oxygen_mode,
        rdf_show_first_shell=args.rdf_show_first_shell,
        rdf_cutoff_mode=args.rdf_cutoff_mode,
        project_mode=args.project_mode,
        coord_cutoff_nm=args.coord_cutoff_nm,
        coord_cutoff_nm_global=args.coord_cutoff_nm_global,
    )
