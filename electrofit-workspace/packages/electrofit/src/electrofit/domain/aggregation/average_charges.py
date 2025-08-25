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
from typing import Tuple, Dict

from electrofit.io.mol2_ops import update_mol2_charges, parse_charges_from_mol2
from electrofit.io.symmetry import load_symmetry_groups
from electrofit.cli.run_commands import run_acpype
from electrofit.config.loader import load_config, dump_config
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

    compose_snapshot(
        results_dir,
        project_root,
        mol_dir.name,
        multi_molecule=multi_mol,
        log_fn=logging.info,
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
    adjust_sym = getattr(proj, "adjust_symmetry", False)
    ignore_sym = getattr(proj, "ignore_symmetry", False)
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
        build_aggregation_decision(
            protocol=protocol,
            adjust_sym=adjust_sym,
            ignore_sym=ignore_sym,
            calc_group_average=calc_group_average,
            group_average_applied=False,
            symmetry_json_found=symmetry_json_found,
        ).log('step6')
        log_relevant_config('step6', proj, ['molecule_name','protocol','adjust_symmetry','ignore_symmetry','calculate_group_average'])
    except Exception:  # pragma: no cover
        logging.debug('[step6][decisions] logging failed', exc_info=True)

    initial_charges_dict = adjust_atom_names(parse_charges_from_mol2(mol2_source_file_path))
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

    try:
        plot_charges_by_atom(atoms_dict, initial_charges_dict, str(results_dir))
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
    cleaned_charges_file: Path | None = None
    executed_acpype = False
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
            current_net = sum(v["average_charge"] for v in cleaned_dict.values())
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
            (results_dir / "cleaned_adjusted_charges.json").write_text(json.dumps(cleaned_dict, indent=2))
            cleaned_lines = [f"{rec['average_charge']:.4f}" for rec in cleaned_dict.values()]
            cleaned_charges_file = results_dir / "cleaned_average_charges.chg"
            cleaned_charges_file.write_text("\n".join(cleaned_lines) + "\n")
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

    # If no group average requested we might still want histograms with just before/after outlier removal.
    # Group averaging logic below may produce adjusted histogram.

    if adjust_sym and not ignore_sym:
        sym_json = extracted / "equiv_groups.json"
        if not sym_json.is_file():
            cand = list(pis_dir.glob("*.json"))
            if cand:
                sym_json = cand[0]
        if sym_json.is_file():
            equiv_group = load_symmetry_groups(str(sym_json))
            try:
                plot_charges_by_symmetry(atoms_dict, initial_charges_dict, str(results_dir), equiv_group)
            except Exception as e:  # pragma: no cover
                logging.debug(f"[step6][symmetry-plot] failed: {e}")
            if calc_group_average and not remove_outlier:
                updated = calculate_symmetric_group_averages(results_dir / "charges_dict.json", sym_json)
                (results_dir / "group_average_charges_dict.json").write_text(json.dumps(updated, indent=2))
                group_txt = results_dir / "group_average_charges.txt"
                with group_txt.open("w") as f:
                    f.write("#Atom_Name\tAverage_Charge\n")
                    for atom, rec in updated.items():
                        f.write(f"{atom}\t{rec['average_charge']:.4f}\n")
                lines = [f"{rec['average_charge']:.4f}" for rec in updated.values()]
                (results_dir / "group_average_charges.chg").write_text("\n".join(lines) + "\n")
                updated_mol2_out = results_dir / f"averaged_{molecule_name}.mol2"
                chg_file = results_dir / "group_average_charges.chg"
                logging.info("[step6] Updating MOL2 with group average charges and running acpype (user charges mode)...")
                update_mol2_charges(mol2_source_file_path, str(chg_file), str(updated_mol2_out))
                run_acpype(str(updated_mol2_out), charge, str(results_dir), proj.atom_type or "gaff2", charges="user")
                try:
                    build_aggregation_decision(
                        protocol=protocol,
                        adjust_sym=adjust_sym,
                        ignore_sym=ignore_sym,
                        calc_group_average=calc_group_average,
                        group_average_applied=True,
                        symmetry_json_found=symmetry_json_found,
                    ).log('step6')
                except Exception:  # pragma: no cover
                    logging.debug('[step6][decisions] post-apply logging failed', exc_info=True)
                if plot_histograms:
                    try:
                        specs_adj = [
                            HistogramSpec(
                                column=a,
                                data_before=atoms_dict[a]["charges"],
                                data_after=None,
                                mean_line=updated[a]["average_charge"],
                                color="darkgreen",
                            )
                            for a in atoms_dict
                        ]
                        plot_atom_histograms(
                            specs_adj,
                            results_dir / "hist_adjusted.pdf",
                            "Charge Distribution with Group Average Charges",
                            bins=hist_bins,
                        )
                        _mark("adjusted", True, True, "hist_adjusted.pdf")
                    except Exception as e:  # pragma: no cover
                        logging.debug(f"[step6][hist] adjusted (group avg) failed: {e}")
                        _mark("adjusted", False, True, None, f"error: {e}")
                if plot_histograms and hist_combine_groups:
                    try:
                        group_specs = []
                        for rep, group_atoms in equiv_group.items():
                            atoms_all = [rep] + list(group_atoms)
                            combined = []
                            for a in atoms_all:
                                combined.extend(atoms_dict.get(a, {}).get("charges", []))
                            mean_line = updated.get(rep, {}).get("average_charge")
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
                        logging.debug(f"[step6][hist] group combined failed: {e}")
                        _mark("group_combined", False, True, None, f"error: {e}")
            else:
                avg_chg = results_dir / "average_charges.chg"
                if avg_chg.is_file() and not calc_group_average and not remove_outlier:  # avoid premature update if group avg path or outlier path will handle it
                    updated_mol2_out = results_dir / f"averaged_{molecule_name}.mol2"
                    logging.info("[step6] Updating MOL2 with average charges (symmetry, no group average) and running acpype...")
                    update_mol2_charges(mol2_source_file_path, str(avg_chg), str(updated_mol2_out))
                    run_acpype(str(updated_mol2_out), charge, str(results_dir), proj.atom_type or "gaff2", charges="user")
                    executed_acpype = True
                elif not avg_chg.is_file():
                    print(f"[step6][warn] (symmetry,no-group-average) missing average_charges.chg in {results_dir}")
                if plot_histograms and hist_combine_groups:
                    # Provide group combined hist even without group average (raw distributions)
                    try:
                        group_specs = []
                        for rep, group_atoms in equiv_group.items():
                            atoms_all = [rep] + list(group_atoms)
                            combined = []
                            for a in atoms_all:
                                combined.extend(atoms_dict.get(a, {}).get("charges", []))
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
                        logging.debug(f"[step6][hist] group combined (no group avg) failed: {e}")
                        _mark("group_combined", False, True, None, f"error: {e}")
        else:
            print(f"[step6][warn] no symmetry JSON for {molecule_name}; skipping symmetry averaging")
            if plot_histograms and hist_combine_groups:
                _mark("group_combined", False, True, None, "missing symmetry JSON")
    else:
        if not calc_group_average and not remove_outlier:  # simple average case, perform mol2 update only when no further adjustments
            avg_chg = results_dir / "average_charges.chg"
            if avg_chg.is_file():
                updated_mol2_out = results_dir / f"averaged_{molecule_name}.mol2"
                logging.info("[step6] Updating MOL2 with average charges (no symmetry averaging) and running acpype...")
                update_mol2_charges(mol2_source_file_path, str(avg_chg), str(updated_mol2_out))
                run_acpype(str(updated_mol2_out), charge, str(results_dir), proj.atom_type or "gaff2", charges="user")
                executed_acpype = True
            else:
                print(f"[step6][warn] expected average_charges.chg missing in {results_dir}")
            if plot_histograms and not remove_outlier:
                # Show average charge means as dashed lines
                try:
                    with open(results_dir / "charges_dict.json") as f:
                        orig = json.load(f)
                    specs_avg = [
                        HistogramSpec(
                            column=a,
                            data_before=orig[a]["charges"],
                            mean_line=orig[a]["average_charge"],
                            color="darkblue",
                        )
                        for a in orig
                    ]
                    plot_atom_histograms(
                        specs_avg,
                        results_dir / "hist_adjusted.pdf",
                        "Charge Distribution with Average Charges",
                        bins=hist_bins,
                    )
                    _mark("adjusted", True, True, "hist_adjusted.pdf")
                except Exception as e:  # pragma: no cover
                    logging.debug(f"[step6][hist] avg lines failed: {e}")
                    _mark("adjusted", False, True, None, f"error: {e}")
    # Outlier + group average combined path (legacy cleaned_adjusted logic)
    # Relaxed condition: if outliers removed and group average requested and a symmetry JSON is available (even if adjust_sym not explicitly set)
    if remove_outlier and calc_group_average:
        # after earlier removal we have cleaned_dict
        try:
            sym_json = extracted / "equiv_groups.json"
            if not sym_json.is_file():
                cand = list((mol_dir / 'run_gau_create_gmx_in').glob('*.json'))
                if cand:
                    sym_json = cand[0]
            if sym_json.is_file() and not ignore_sym:
                updated = calculate_symmetric_group_averages(results_dir / "cleaned_adjusted_charges.json", sym_json)
                (results_dir / "cleaned_adjusted_group_average_charges_dict.json").write_text(json.dumps(updated, indent=2))
                lines = [f"{rec['average_charge']:.4f}" for rec in updated.values()]
                cleaned_group_chg = results_dir / "cleaned_adjusted_group_average_charges.chg"
                cleaned_group_chg.write_text("\n".join(lines) + "\n")
                # Update MOL2 + run acpype immediately for this combined path
                try:
                    updated_mol2_out = results_dir / f"averaged_{molecule_name}_cleaned.mol2"
                    update_mol2_charges(mol2_source_file_path, str(cleaned_group_chg), str(updated_mol2_out))
                    run_acpype(str(updated_mol2_out), charge, str(results_dir), proj.atom_type or "gaff2", charges="user")
                    executed_acpype = True
                except Exception:  # pragma: no cover
                    logging.debug("[step6][outlier+group] mol2/acpype failed", exc_info=True)
                if plot_histograms:
                    try:
                        specs_clean_adj = [
                            HistogramSpec(
                                column=a,
                                data_before=atoms_dict[a]["charges"],
                                data_after=cleaned_dict[a]["charges"],
                                mean_line=updated[a]["average_charge"],
                                color="darkred",
                                overlay_color="darkgreen",
                            )
                            for a in atoms_dict
                        ]
                        plot_atom_histograms(
                            specs_clean_adj,
                            results_dir / "hist_adjusted.pdf",
                            "Charge Distribution of Clipped Data with Reweighted Group Average Charges",
                            bins=hist_bins,
                        )
                        _mark("adjusted", True, True, "hist_adjusted.pdf")
                    except Exception as e:  # pragma: no cover
                        logging.debug(f"[step6][hist] cleaned adjusted failed: {e}")
                        _mark("adjusted", False, True, None, f"error: {e}")
                    if hist_combine_groups:
                        try:
                            equiv_group = load_symmetry_groups(str(sym_json))
                            group_specs = []
                            for rep, group_atoms in equiv_group.items():
                                atoms_all = [rep] + list(group_atoms)
                                combined_before = []
                                combined_after = []
                                for a in atoms_all:
                                    combined_before.extend(atoms_dict.get(a, {}).get("charges", []))
                                    combined_after.extend(cleaned_dict.get(a, {}).get("charges", []))
                                mean_line = updated.get(rep, {}).get("average_charge") if isinstance(updated.get(rep, {}), dict) else None
                                group_specs.append(
                                    HistogramSpec(
                                        column=f"Group:{rep}",
                                        data_before=combined_before,
                                        data_after=combined_after,
                                        mean_line=mean_line,
                                        color="purple",
                                        overlay_color="darkgreen",
                                    )
                                )
                            if group_specs:
                                plot_atom_histograms(
                                    group_specs,
                                    results_dir / "hist_groups.pdf",
                                    "Group-Combined Charge Distributions (Cleaned)",
                                    bins=hist_bins,
                                )
                                _mark("group_combined", True, True, "hist_groups.pdf")
                            else:
                                _mark("group_combined", False, True, None, "no groups")
                        except Exception as e:  # pragma: no cover
                            logging.debug(f"[step6][hist] cleaned group combined failed: {e}")
                            _mark("group_combined", False, True, None, f"error: {e}")
        except Exception as e:  # pragma: no cover
            logging.debug(f"[step6][outlier+group] failed: {e}")
            if plot_histograms:
                _mark("adjusted", False, True, None, f"error: {e}")
                if hist_combine_groups:
                    _mark("group_combined", False, True, None, f"error: {e}")
    # Post-outlier finalisation: if outlier_removed and acpype noch nicht lief, nutze gereinigte oder group-average Dateien.
    if outlier_removed and not executed_acpype:
        candidate = None
        out_name = None
        if (results_dir / "cleaned_adjusted_group_average_charges.chg").is_file():
            candidate = results_dir / "cleaned_adjusted_group_average_charges.chg"
            out_name = f"averaged_{molecule_name}_cleaned.mol2"
        elif cleaned_charges_file and cleaned_charges_file.is_file():
            candidate = cleaned_charges_file
            out_name = f"averaged_{molecule_name}_cleaned.mol2"
        elif (results_dir / "group_average_charges.chg").is_file():
            candidate = results_dir / "group_average_charges.chg"
            out_name = f"averaged_{molecule_name}.mol2"
        elif (results_dir / "average_charges.chg").is_file():
            candidate = results_dir / "average_charges.chg"
            out_name = f"averaged_{molecule_name}.mol2"
        if candidate and out_name:
            try:
                updated_mol2_out = results_dir / out_name
                logging.info("[step6] Finalizing MOL2 update after outlier removal and running acpype...")
                update_mol2_charges(mol2_source_file_path, str(candidate), str(updated_mol2_out))
                run_acpype(str(updated_mol2_out), charge, str(results_dir), proj.atom_type or "gaff2", charges="user")
            except Exception:  # pragma: no cover
                logging.debug("[step6] post-outlier mol2/acpype failed", exc_info=True)
    # Write histogram manifest
    if plot_histograms:
        try:
            (results_dir / "hist_manifest.json").write_text(json.dumps(hist_manifest, indent=2))
        except Exception:  # pragma: no cover
            logging.debug("[step6][hist] failed writing manifest", exc_info=True)
    return True, "ok"
