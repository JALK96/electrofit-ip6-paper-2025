"""Domain logic for Step6: aggregation of RESP ensemble charges and symmetry/group averaging.

Separated from the pipeline step so that orchestration code stays thin and
testable; this module contains only deterministic transformation logic plus
plotting side–effects.
"""
from __future__ import annotations

from pathlib import Path
import fnmatch
import json
import logging
import os
import math
import statistics
from typing import Tuple, Dict

from electrofit.io.mol2_ops import update_mol2_charges, parse_charges_from_mol2
from electrofit.io.symmetry import load_symmetry_groups
from electrofit.cli.run_commands import run_acpype
from electrofit.config.loader import load_config, dump_config, resolve_symmetry_flags
from electrofit.infra.config_snapshot import compose_snapshot
from electrofit.io.files import (
    adjust_atom_names,
    extract_charges_from_subdirectories,
)

from electrofit.viz.helpers import plot_charges_by_symmetry, plot_charges_by_atom
from electrofit.viz.histograms import HistogramSpec, plot_atom_histograms
from electrofit.infra.step_logging import log_relevant_config
from electrofit.infra.decisions import build_aggregation_decision

__all__ = ["process_molecule_average_charges", "calculate_symmetric_group_averages"]


def calculate_symmetric_group_averages(charges_dict_file: Path, equivalent_groups_file: Path):
    """Return updated charges dict where all atoms in a symmetry group share the mean charge."""
    with open(equivalent_groups_file, "r") as f:
        equivalent_groups = json.load(f)
    with open(charges_dict_file, "r") as f:
        charges_dict = json.load(f)
    updated = charges_dict.copy()
    for representative, group in equivalent_groups.items():
        full_group = [representative] + list(group)
        vals = [charges_dict[a]["average_charge"] for a in full_group if a in charges_dict]
        if not vals:
            continue
        mean_v = sum(vals) / len(vals)
        for atom in full_group:
            if atom in updated:
                updated[atom]["average_charge"] = mean_v
    return updated


def _apply_group_average_in_memory(charges_dict: dict[str, dict], equivalent_groups: dict[str, list[str]]) -> dict[str, dict]:
    """Return a deep-copied charges dict with symmetry-group averages applied."""
    updated = {atom: {"charges": list(rec.get("charges", [])), "average_charge": rec.get("average_charge", 0.0)} for atom, rec in charges_dict.items()}
    for representative, group in equivalent_groups.items():
        full_group = [representative] + list(group)
        vals = [charges_dict.get(a, {}).get("average_charge") for a in full_group if a in charges_dict]
        vals = [v for v in vals if v is not None]
        if not vals:
            continue
        mean_v = sum(vals) / len(vals)
        for atom in full_group:
            if atom in updated:
                updated[atom]["average_charge"] = mean_v
    return updated


def _write_charge_bundle(results_dir: Path, prefix: str, charges_dict: dict[str, dict]) -> Path:
    """Persist charges dict in .json/.txt/.chg formats and return the .chg path."""
    json_path = results_dir / f"{prefix}_dict.json"
    json_path.write_text(json.dumps(charges_dict, indent=2))
    txt_path = results_dir / f"{prefix}.txt"
    with txt_path.open("w") as f:
        f.write("#Atom_Name\tAverage_Charge\n")
        for atom, rec in charges_dict.items():
            f.write(f"{atom}\t{rec['average_charge']:.4f}\n")
    chg_path = results_dir / f"{prefix}.chg"
    chg_lines = [f"{rec['average_charge']:.4f}" for rec in charges_dict.values()]
    chg_path.write_text("\n".join(chg_lines) + "\n")
    return chg_path


def _find_symmetry_json(results_dir: Path, extracted: Path, pis_dir: Path) -> Path | None:
    """Locate a symmetry JSON file, checking results, extracted_conforms and run directories."""
    candidates: list[Path] = [
        results_dir / "equiv_groups.json",
        extracted / "equiv_groups.json",
    ]
    candidates.extend(sorted(pis_dir.glob("*.json")))
    for cand in candidates:
        if cand.is_file():
            return cand
    return None


def process_molecule_average_charges(
    mol_dir: Path,
    project_root: Path,
    override_cfg: Path | None,
    multi_mol: bool,
    remove_outlier: bool,
    plot_histograms: bool = False,
    hist_combine_groups: bool = False,
    hist_bins: int = 20,
    outlier_iqr_factor: float = 1.5,
) -> Tuple[bool, str]:
    """Aggregate charges for a single molecule directory.

    Returns (ok, message) mirroring legacy behaviour for backwards compatible console output.

    Notes
    -----
    Experimental: The parameter ``remove_outlier`` enables an IQR based conformer
    filtering pass. This is still experimental: (a) the statistical robustness
    for small conformer ensembles (< ~15) is weak, (b) downstream re‑normalisation
    of net charge can amplify noise on sparsely sampled atoms, and (c) symmetry /
    group averaging after heavy pruning may bias distributions. Use only for
    exploratory analysis; not yet recommended for production charge sets.
    """
    extracted = mol_dir / "extracted_conforms"
    if remove_outlier:
        logging.warning(
            "[step6][experimental] Outlier removal is EXPERIMENTAL and not yet recommended to use."
        )
    pis_dir = mol_dir / "run_gau_create_gmx_in"
    results_dir = mol_dir / "results"
    if not extracted.is_dir():
        return False, "no extracted_conforms"
    acpype_glob = list(pis_dir.glob("*.acpype"))
    if not acpype_glob:
        return False, "no acpype dir"
    ac_dir = acpype_glob[0]
    results_dir.mkdir(exist_ok=True)
    # Plot artifacts go into a dedicated subfolder to avoid clobbering root-level *.chg outputs
    plots_dir = results_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    compose_snapshot(
        results_dir,
        project_root,
        mol_dir.name,
        multi_molecule=multi_mol,
        log_fn=logging.info,
        step="step6",
        upstream=extracted / "electrofit.toml",
        process_cfg=mol_dir / "electrofit.toml",
        molecule_input=project_root / "data" / "input" / mol_dir.name / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
        extra_override=override_cfg,
    )
    cfg = load_config(project_root, context_dir=results_dir, molecule_name=mol_dir.name)
    try:
        dump_config(cfg, log_fn=logging.debug)
    except Exception as e:  # pragma: no cover
        logging.debug(f"[step6] config dump failed: {e}")

    proj = cfg.project
    molecule_name = proj.molecule_name or mol_dir.name
    charge = proj.charge or 0
    atom_type = proj.atom_type or "gaff2"
    adjust_sym, ignore_sym = resolve_symmetry_flags(cfg, "ensemble")
    proj.adjust_symmetry = adjust_sym  # type: ignore[attr-defined]
    proj.ignore_symmetry = ignore_sym  # type: ignore[attr-defined]
    calc_group_average = getattr(proj, "calculate_group_average", False)
    protocol = getattr(proj, "protocol", None)

    mol2_source_file_path = None
    pattern = f"*{atom_type}.mol2"
    for fn in os.listdir(ac_dir):
        if fnmatch.fnmatch(fn, pattern):
            mol2_source_file_path = os.path.join(ac_dir, fn)
            break
    if not mol2_source_file_path:
        return False, "no matching mol2"

    try:
        symmetry_json_found = (extracted / 'equiv_groups.json').is_file() or any((mol_dir / 'run_gau_create_gmx_in').glob('*.json'))
        log_relevant_config(
            'step6',
            cfg,
            [
                'project.molecule_name',
                'project.protocol',
                'symmetry.initial',
                'symmetry.ensemble',
                'project.adjust_symmetry',
                'project.ignore_symmetry',
                'project.calculate_group_average',
            ],
        )
    except Exception:  # pragma: no cover
        logging.debug('[step6][decisions] logging failed', exc_info=True)

    initial_charges_dict = adjust_atom_names(parse_charges_from_mol2(mol2_source_file_path))

    # ── Step6 input audit: count conformers + sanity-check per-conformer RESP outputs ─────────
    # This helps detect missing/failed conformers (no *_resp.mol2) and pathological outputs
    # (e.g., all-zero or non-finite charges) early and transparently.
    try:
        target_net = float(charge)
        conf_dirs = sorted([p for p in extracted.iterdir() if p.is_dir()])
        missing_resp: list[str] = []
        suspicious_resp: list[str] = []
        resp_nets: list[float] = []
        resp_atom_counts: list[int] = []
        large_net_deviation: list[str] = []
        parse_fail: list[str] = []

        for conf_dir in conf_dirs:
            resp_files = sorted(conf_dir.glob("*_resp.mol2"))
            if not resp_files:
                missing_resp.append(conf_dir.name)
                continue
            resp_path = resp_files[0]
            try:
                atoms = parse_charges_from_mol2(str(resp_path))
                charges = []
                for rec in atoms.values():
                    vals = rec.get("charges") or []
                    if vals:
                        charges.append(float(vals[0]))
                if not charges:
                    suspicious_resp.append(conf_dir.name)
                    continue
                if any(not math.isfinite(v) for v in charges):
                    suspicious_resp.append(conf_dir.name)
                    continue
                if max(abs(v) for v in charges) < 1e-8:
                    # Classic failure mode: everything is 0.0 in *_resp.mol2
                    suspicious_resp.append(conf_dir.name)
                net = float(sum(charges))
                resp_nets.append(net)
                resp_atom_counts.append(len(charges))
                if abs(net - target_net) > 0.5:
                    large_net_deviation.append(conf_dir.name)
            except Exception:
                parse_fail.append(conf_dir.name)

        total_dirs = len(conf_dirs)
        with_resp = total_dirs - len(missing_resp)
        logging.info(
            "[step6] Conformer RESP outputs: dirs=%d with_resp=%d missing_resp=%d suspicious=%d parse_fail=%d",
            total_dirs,
            with_resp,
            len(missing_resp),
            len(suspicious_resp),
            len(parse_fail),
        )
        if missing_resp:
            logging.info("[step6] Missing *_resp.mol2 (first 10): %s", ", ".join(missing_resp[:10]))
        if suspicious_resp:
            logging.warning("[step6] Suspicious *_resp.mol2 (first 10): %s", ", ".join(suspicious_resp[:10]))
        if parse_fail:
            logging.warning("[step6] Failed to parse *_resp.mol2 (first 10): %s", ", ".join(parse_fail[:10]))
        if resp_nets:
            mean_net = statistics.fmean(resp_nets)
            std_net = statistics.pstdev(resp_nets) if len(resp_nets) > 1 else 0.0
            min_net = min(resp_nets)
            max_net = max(resp_nets)
            logging.info(
                "[step6] Per-conformer net charge: mean=%.4f std=%.4f min=%.4f max=%.4f target=%.1f",
                mean_net,
                std_net,
                min_net,
                max_net,
                target_net,
            )
        if large_net_deviation:
            logging.warning(
                "[step6] Conformers with |net-target| > 0.5 e (first 10): %s",
                ", ".join(large_net_deviation[:10]),
            )
        if resp_atom_counts:
            exp_atoms = len(initial_charges_dict)
            # If atom counts vary, this usually indicates corrupted/mol2 parsing issues.
            if any(n != exp_atoms for n in resp_atom_counts):
                bad = sum(1 for n in resp_atom_counts if n != exp_atoms)
                logging.warning(
                    "[step6] Atom count mismatch in *_resp.mol2: expected=%d; mismatching_conformers=%d",
                    exp_atoms,
                    bad,
                )
    except Exception:  # pragma: no cover
        logging.debug("[step6] input audit failed", exc_info=True)

    atoms_dict = extract_charges_from_subdirectories(str(extracted), str(results_dir))
    (results_dir / "charges_dict.json").write_text(json.dumps(atoms_dict, indent=2))
    (results_dir / "initial_charges_dict.json").write_text(json.dumps(initial_charges_dict, indent=2))

    # Ensure we always emit average_charges.chg early (robustness for downstream mol2 update)
    try:
        avg_lines = []
        for atom, rec in atoms_dict.items():
            ch_list = rec.get("charges", []) or []
            avg = (sum(ch_list) / len(ch_list)) if ch_list else 0.0
            rec["average_charge"] = avg
            avg_lines.append(f"{avg:.4f}")
        if avg_lines:
            (results_dir / "average_charges.chg").write_text("\n".join(avg_lines) + "\n")
        else:
            # fallback: replicate from average_charges.txt if that exists later
            logging.warning("[step6] no charge data collected; skipping average_charges.chg emission")
    except Exception:  # pragma: no cover
        logging.debug("[step6] early average_charges.chg write failed", exc_info=True)

    # Net charge summary for the raw mean charges (before any optional filtering/group averaging).
    try:
        target_net = float(charge)
        raw_net = float(sum(rec.get("average_charge", 0.0) for rec in atoms_dict.values()))
        delta = raw_net - target_net
        logging.info(
            "[step6] Net charge (raw mean charges): net=%.4f target=%.1f Δ=%.4f",
            raw_net,
            target_net,
            delta,
        )
        if abs(delta) > 0.1:
            logging.warning("[step6] Net charge deviates from target by > 0.1 e (raw mean charges).")
    except Exception:  # pragma: no cover
        logging.debug("[step6] net charge (raw) summary failed", exc_info=True)

    try:
        plot_charges_by_atom(atoms_dict, initial_charges_dict, str(plots_dir))
    except Exception as e:  # pragma: no cover
        print(f"[step6][warn] plotting per-atom charges failed: {e}")

    # Histogram manifest collects status of each attempted histogram
    hist_manifest = {}
    def _mark(name: str, created: bool, expected: bool, path: str | None, reason: str | None = None):
        hist_manifest[name] = {
            "expected": expected,
            "created": created,
            "path": path,
            "reason": reason,
        }

    # Optional histogram BEFORE any filtering
    if plot_histograms:
        try:
            specs = []
            for atom, rec in atoms_dict.items():
                data = [v for v in rec.get("charges", []) if v is not None]
                if not data:
                    continue  # skip empty distributions
                specs.append(
                    HistogramSpec(
                        column=atom,
                        data_before=data,
                        color="darkred",
                    )
                )
            if specs:
                target = results_dir / "hist.pdf"
                plot_atom_histograms(specs, target, "Charge Distribution Before Filtering", bins=hist_bins)
                created = target.is_file()
                _mark("initial", created, True, "hist.pdf" if created else None, None if created else "not created")
            else:
                logging.info("[step6][hist] no non-empty charge series for initial histogram")
                _mark("initial", False, True, None, "no non-empty series")
        except Exception as e:  # pragma: no cover
            logging.debug(f"[step6][hist] failed initial hist plot: {e}")
            _mark("initial", False, True, None, f"error: {e}")
    else:
        _mark("initial", False, False, None, "plot_histograms flag not set")

    cleaned_dict = atoms_dict
    outlier_removed = False
    if remove_outlier:
        # Construct dataframe for IQR filtering (row = conformer)
        try:
            import pandas as pd  # local import to avoid hard dep at import time
            df = pd.DataFrame({k: v["charges"] for k, v in atoms_dict.items()})
            # build per-row mask: mark rows where ANY atom charge is outside IQR*factor
            def iqr_mask(series, factor):
                q1 = series.quantile(0.25)
                q3 = series.quantile(0.75)
                iqr = q3 - q1
                return (series < (q1 - factor * iqr)) | (series > (q3 + factor * iqr))
            row_outlier = pd.Series(False, index=df.index)
            per_atom_outliers: Dict[str, int] = {}
            for col in df.columns:
                col_mask = iqr_mask(df[col], outlier_iqr_factor)
                per_atom_outliers[col] = int(col_mask.sum())
                row_outlier = row_outlier | col_mask
            df_clean = df[~row_outlier].reset_index(drop=True)
            if plot_histograms:
                specs_no = [
                    HistogramSpec(
                        column=col,
                        data_before=df[col].dropna().tolist(),
                        data_after=df_clean[col].dropna().tolist(),
                        color="darkred",
                        overlay_color="darkblue",
                    )
                    for col in df.columns
                ]
                target2 = results_dir / "hist_no_outlier.pdf"
                plot_atom_histograms(
                    specs_no,
                    target2,
                    "Charge Distribution After Outlier Removal",
                    bins=hist_bins,
                )
                created2 = target2.is_file()
                _mark("after_outlier", created2, True, "hist_no_outlier.pdf" if created2 else None, None if created2 else "not created")
            # reconstruct cleaned_dict structure
            cleaned_dict = {}
            for atom in df_clean.columns:
                charges_list = df_clean[atom].dropna().tolist()
                avg = sum(charges_list) / len(charges_list) if charges_list else 0.0
                cleaned_dict[atom] = {"charges": charges_list, "average_charge": avg}
            # Net charge normalisation (legacy parity + improvement)
            try:
                target_net = float(charge)
            except Exception:
                target_net = 0.0
            current_net = float(sum(v["average_charge"] for v in cleaned_dict.values()))
            deviation = target_net - current_net
            if cleaned_dict:
                if abs(current_net) < 1e-12:
                    per_atom = deviation / len(cleaned_dict)
                    for rec in cleaned_dict.values():
                        rec["average_charge"] += per_atom
                else:
                    for rec in cleaned_dict.values():
                        proportion = rec["average_charge"] / current_net if current_net != 0 else 0.0
                        rec["average_charge"] += proportion * deviation
            try:
                after_net = float(sum(v["average_charge"] for v in cleaned_dict.values()))
                logging.info(
                    "[step6] Net charge (cleaned): before=%.4f after=%.4f target=%.1f",
                    current_net,
                    after_net,
                    target_net,
                )
            except Exception:
                pass
            (results_dir / "cleaned_adjusted_charges.json").write_text(json.dumps(cleaned_dict, indent=2))
            cleaned_lines = [f"{rec['average_charge']:.4f}" for rec in cleaned_dict.values()]
            cleaned_charges_path = results_dir / "cleaned_average_charges.chg"
            cleaned_charges_path.write_text("\n".join(cleaned_lines) + "\n")
            outlier_removed = True
            logging.info(
                f"[step6] Outlier removal: removed {row_outlier.sum()} / {len(row_outlier)} conformers (factor={outlier_iqr_factor})."
            )
            summary = {
                "total_conformers": int(len(row_outlier)),
                "removed_conformers": int(row_outlier.sum()),
                "kept_conformers": int(len(row_outlier) - row_outlier.sum()),
                "iqr_factor": outlier_iqr_factor,
                "per_atom_outliers": per_atom_outliers,
                "removed_indices": [int(i) for i, v in row_outlier.items() if v],
            }
            (results_dir / "hist_summary.json").write_text(json.dumps(summary, indent=2))
        except Exception as e:  # pragma: no cover
            logging.debug(f"[step6][outlier] removal failed: {e}")
            if plot_histograms:
                _mark("after_outlier", False, True, "hist_no_outlier.pdf", f"error: {e}")

    # Determine symmetry mapping and optional group averaging behaviour
    sym_json_path = _find_symmetry_json(results_dir, extracted, pis_dir)
    symmetry_mapping: dict[str, list[str]] | None = None
    symmetry_json_found = False
    if sym_json_path:
        try:
            symmetry_mapping = load_symmetry_groups(str(sym_json_path))
            symmetry_json_found = True
        except Exception as e:  # pragma: no cover
            logging.debug(f"[step6][symmetry] failed to load {sym_json_path}: {e}")
    elif calc_group_average:
        logging.warning("[step6] calculate_group_average requested but no symmetry JSON found; skipping group averaging")

    base_dict_for_group = cleaned_dict if outlier_removed else atoms_dict

    if symmetry_mapping:
        try:
            plot_charges_by_symmetry(atoms_dict, initial_charges_dict, str(plots_dir), symmetry_mapping)
        except Exception as e:  # pragma: no cover
            logging.debug(f"[step6][symmetry-plot] failed: {e}")
        try:
            plot_charges_by_atom_sym(atoms_dict, initial_charges_dict, str(plots_dir), equivalent_groups=symmetry_mapping)
        except Exception as e:  # pragma: no cover
            logging.debug(f"[step6][symmetry-plot-atom] failed: {e}")
        if plot_histograms and hist_combine_groups and not calc_group_average:
            try:
                group_specs = []
                for rep, group_atoms in symmetry_mapping.items():
                    atoms_all = [rep] + list(group_atoms)
                    combined = []
                    for atom in atoms_all:
                        combined.extend(base_dict_for_group.get(atom, {}).get("charges", []))
                    mean_line = sum(combined) / len(combined) if combined else None
                    group_specs.append(
                        HistogramSpec(
                            column=f"Group:{rep}",
                            data_before=combined,
                            mean_line=mean_line,
                            color="purple",
                        )
                    )
                if group_specs:
                    plot_atom_histograms(
                        group_specs,
                        results_dir / "hist_groups.pdf",
                        "Group-Combined Charge Distributions",
                        bins=hist_bins,
                    )
                    _mark("group_combined", True, True, "hist_groups.pdf")
                else:
                    _mark("group_combined", False, True, None, "no groups")
            except Exception as e:  # pragma: no cover
                logging.debug(f"[step6][hist] group combined (raw) failed: {e}")
                _mark("group_combined", False, True, None, f"error: {e}")
    else:
        if calc_group_average:
            _mark("group_combined", False, True, None, "missing symmetry JSON")

    group_avg_dict: dict[str, dict] | None = None
    group_avg_chg_path: Path | None = None
    group_average_applied = False
    if calc_group_average and symmetry_mapping:
        group_avg_dict = _apply_group_average_in_memory(base_dict_for_group, symmetry_mapping)
        prefix = "cleaned_adjusted_group_average_charges" if outlier_removed else "group_average_charges"
        group_avg_chg_path = _write_charge_bundle(results_dir, prefix, group_avg_dict)
        group_average_applied = True
        try:
            target_net = float(charge)
            group_net = float(sum(rec.get("average_charge", 0.0) for rec in group_avg_dict.values()))
            logging.info(
                "[step6] Net charge (group-averaged): net=%.4f target=%.1f Δ=%.4f",
                group_net,
                target_net,
                group_net - target_net,
            )
        except Exception:  # pragma: no cover
            logging.debug("[step6] net charge (group avg) summary failed", exc_info=True)
        if plot_histograms:
            try:
                specs_adj = [
                    HistogramSpec(
                        column=atom,
                        data_before=base_dict_for_group[atom]["charges"],
                        data_after=None,
                        mean_line=group_avg_dict[atom]["average_charge"],
                        color="darkgreen",
                    )
                    for atom in base_dict_for_group
                ]
                plot_atom_histograms(
                    specs_adj,
                    results_dir / "hist_adjusted.pdf",
                    "Charge Distribution with Group Average Charges",
                    bins=hist_bins,
                )
                _mark("adjusted", True, True, "hist_adjusted.pdf")
            except Exception as e:  # pragma: no cover
                logging.debug(f"[step6][hist] group average means failed: {e}")
                _mark("adjusted", False, True, None, f"error: {e}")
        if plot_histograms and hist_combine_groups:
            try:
                group_specs = []
                for rep, group_atoms in symmetry_mapping.items():
                    atoms_all = [rep] + list(group_atoms)
                    combined = []
                    for atom in atoms_all:
                        combined.extend(base_dict_for_group.get(atom, {}).get("charges", []))
                    mean_line = group_avg_dict.get(rep, {}).get("average_charge")
                    group_specs.append(
                        HistogramSpec(
                            column=f"Group:{rep}",
                            data_before=combined,
                            mean_line=mean_line,
                            color="purple",
                        )
                    )
                if group_specs:
                    plot_atom_histograms(
                        group_specs,
                        results_dir / "hist_groups.pdf",
                        "Group-Combined Charge Distributions",
                        bins=hist_bins,
                    )
                    _mark("group_combined", True, True, "hist_groups.pdf")
                else:
                    _mark("group_combined", False, True, None, "no groups")
            except Exception as e:  # pragma: no cover
                logging.debug(f"[step6][hist] group combined (avg) failed: {e}")
                _mark("group_combined", False, True, None, f"error: {e}")
        try:
            plot_charges_by_symmetry(group_avg_dict, initial_charges_dict, str(plots_dir), symmetry_mapping)
        except Exception as e:  # pragma: no cover
            logging.debug(f"[step6][symmetry-plot][avg] failed: {e}")
        try:
            plot_charges_by_atom_sym(
                base_dict_for_group,
                initial_charges_dict,
                str(plots_dir),
                atoms_dict2=group_avg_dict,
                equivalent_groups=symmetry_mapping,
            )
        except Exception as e:  # pragma: no cover
            logging.debug(f"[step6][symmetry-plot-atom][avg] failed: {e}")
    elif plot_histograms and not remove_outlier:
        # Show ensemble averages as reference lines when no group averaging is requested
        try:
            specs_avg = [
                HistogramSpec(
                    column=atom,
                    data_before=atoms_dict[atom]["charges"],
                    mean_line=atoms_dict[atom]["average_charge"],
                    color="darkblue",
                )
                for atom in atoms_dict
            ]
            plot_atom_histograms(
                specs_avg,
                results_dir / "hist_adjusted.pdf",
                "Charge Distribution with Average Charges",
                bins=hist_bins,
            )
            _mark("adjusted", True, True, "hist_adjusted.pdf")
        except Exception as e:  # pragma: no cover
            logging.debug(f"[step6][hist] avg reference failed: {e}")
            _mark("adjusted", False, True, None, f"error: {e}")

    try:
        sym_ensemble = getattr(getattr(cfg, "symmetry", None), "ensemble", None)
        dec = build_aggregation_decision(
            protocol=protocol,
            adjust_sym=adjust_sym,
            ignore_sym=ignore_sym,
            calc_group_average=calc_group_average,
            group_average_applied=group_average_applied,
            symmetry_json_found=symmetry_json_found,
            symmetry_mode=sym_ensemble,
        )
        dec.log('step6')
    except Exception:  # pragma: no cover
        logging.debug('[step6][decisions] aggregation logging failed', exc_info=True)

    # Select final charge file for MOL2 update / ACPYPE
    def _sum_chg(path: Path) -> tuple[float, int]:
        total = 0.0
        n = 0
        for raw in path.read_text().splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            total += float(line)
            n += 1
        return total, n

    charge_candidates: list[tuple[Path, str]] = [
        (results_dir / "cleaned_adjusted_group_average_charges.chg", f"averaged_{molecule_name}_cleaned.mol2"),
        (results_dir / "cleaned_average_charges.chg", f"averaged_{molecule_name}_cleaned.mol2"),
        (results_dir / "group_average_charges.chg", f"averaged_{molecule_name}.mol2"),
        (results_dir / "average_charges.chg", f"averaged_{molecule_name}.mol2"),
    ]
    for chg_path, mol2_name in charge_candidates:
        if chg_path.is_file():
            try:
                target_net = float(charge)
                net, nvals = _sum_chg(chg_path)
                exp_atoms = len(initial_charges_dict)
                logging.info(
                    "[step6] Selected charge set: %s (n=%d expected=%d net=%.4f target=%.1f Δ=%.4f)",
                    chg_path.name,
                    nvals,
                    exp_atoms,
                    net,
                    target_net,
                    net - target_net,
                )
                if nvals != exp_atoms:
                    logging.warning(
                        "[step6] Charge file length mismatch: %s has n=%d, expected=%d (from mol2).",
                        chg_path.name,
                        nvals,
                        exp_atoms,
                    )
                updated_mol2_out = results_dir / mol2_name
                logging.info("[step6] Updating MOL2 with %s and running acpype...", chg_path.name)
                update_mol2_charges(mol2_source_file_path, str(chg_path), str(updated_mol2_out))
                run_acpype(str(updated_mol2_out), charge, str(results_dir), proj.atom_type or "gaff2", charges="user")
            except Exception:  # pragma: no cover
                logging.debug("[step6] mol2/acpype update failed", exc_info=True)
            break

    # Write histogram manifest
    if plot_histograms:
        try:
            (results_dir / "hist_manifest.json").write_text(json.dumps(hist_manifest, indent=2))
        except Exception:  # pragma: no cover
            logging.debug("[step6][hist] failed writing manifest", exc_info=True)
    return True, "ok"
