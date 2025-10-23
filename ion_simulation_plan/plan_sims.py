from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from conc2time import (  # type: ignore
    AVOGADRO,
    _density_stats,
    _filter_by_ion,
    _load_microstate_config,
    _load_records,
    _normalize_ion_symbol,
    _performance_stats,
    _resolve_ion_charge,
    _derive_box_lengths,
)


# Water models and relative runtime multipliers (TIP4P/2005 assumed 33% slower)
WATER_MODELS: Sequence[tuple[str, float]] = (
    ("TIP3P", 1.0),
    ("TIP4P/2005", 4.0 / 3.0),
)

# Ion planning matrix: concentrations given in millimolar
ION_MATRIX = {
    "Na+": {
        "symbol": "NA",
        "concentrations_mM": (10.0, 150.0, 300.0),
    },
    "Ca2+": {
        "symbol": "CA",
        "concentrations_mM": (5.0, 10.0, 20.0),
    },
    "Mg2+": {
        "symbol": "MG",
        "concentrations_mM": (5.0, 10.0, 20.0),
    },
}


ION_MATRIX = {
    "Mg2+": {
        "symbol": "MG",
        "concentrations_mM": (50, 100, 150.0),
    },
}



@dataclass
class RunEstimate:
    ion_label: str
    ion_symbol: str
    concentration_mM: float
    positive_count: int
    negative_count: int
    water_model: str
    multiplier: float
    wall_hours: float
    wall_hours_slow: float
    wall_hours_fast: float
    atoms: int
    box_edge_nm: float
    fallback_used: bool


def _format_number(value: float) -> str:
    if value >= 1e6:
        return f"{value/1e6:.2f}M"
    if value >= 1e3:
        return f"{value/1e3:.1f}k"
    return f"{value:.0f}"


def _format_hours(hours: float) -> str:
    days = hours / 24.0
    if hours < 24.0:
        return f"{hours:.1f} h"
    return f"{hours:.1f} h (~{days:.1f} d)"


def _ensure_neutral_counts(
    positive_count: int,
    cation_charge: int,
    anion_charge: int,
    solute_charge: int,
) -> tuple[int, int]:
    if cation_charge <= 0:
        raise ValueError("Cation charge must be positive.")
    if anion_charge >= 0:
        raise ValueError("Anion charge must be negative.")

    min_positive = math.ceil(-solute_charge / cation_charge) if cation_charge > 0 else 0
    if positive_count < min_positive:
        positive_count = min_positive

    modulus = abs(anion_charge)
    neutrality_term = positive_count * cation_charge + solute_charge
    remainder = neutrality_term % modulus
    if remainder != 0:
        for delta in range(1, modulus + 1):
            test_term = (positive_count + delta) * cation_charge + solute_charge
            if test_term % modulus == 0:
                positive_count += delta
                neutrality_term = test_term
                break
        else:  # pragma: no cover - defensive; should not happen for monovalent anions
            raise ValueError("Unable to satisfy neutrality with integer ion counts.")

    total_positive_charge = positive_count * cation_charge
    if total_positive_charge + solute_charge < 0:
        raise ValueError("Insufficient positive ions to neutralize the solute charge.")

    negative_count = (total_positive_charge + solute_charge) // abs(anion_charge)
    return positive_count, negative_count


def _compute_base_estimate(
    ion_symbol: str,
    concentration_mM: float,
    *,
    sim_ns: float,
    microstate_config,
    density_stats,
    perf_stats,
    cation_charge: int,
    anion_symbol: str,
    anion_charge: int,
) -> tuple[float, float, float, int, float]:
    concentration_mol = concentration_mM / 1000.0
    if concentration_mol <= 0.0:
        raise ValueError("Concentration must be positive.")

    positive_count = microstate_config.positive_count
    if microstate_config.enforce_neutrality:
        positive_count, negative_count = _ensure_neutral_counts(
            positive_count,
            cation_charge,
            anion_charge,
            microstate_config.charge,
        )
    else:
        negative_count = microstate_config.negative_count or 0

    mol_positive = positive_count / AVOGADRO
    volume_litre = mol_positive / concentration_mol
    volume_nm3 = volume_litre * 1e24
    lengths_nm = _derive_box_lengths(volume_nm3, microstate_config.box_type)

    atoms_median = density_stats.median * volume_nm3
    atoms_q1 = density_stats.q1 * volume_nm3
    atoms_q3 = density_stats.q3 * volume_nm3

    ns_day_median = perf_stats.median_k / atoms_median
    ns_day_slow = perf_stats.q1_k / atoms_q3 if atoms_q3 > 0 else float("nan")
    ns_day_fast = perf_stats.q3_k / atoms_q1 if atoms_q1 > 0 else float("nan")

    wall_days_median = sim_ns / ns_day_median
    wall_days_slow = sim_ns / ns_day_slow if ns_day_slow > 0 else float("nan")
    wall_days_fast = sim_ns / ns_day_fast if ns_day_fast > 0 else float("nan")

    box_edge = lengths_nm[0] if lengths_nm else float("nan")
    atoms_estimate = int(round(atoms_median))

    return (
        wall_days_median * 24.0,
        wall_days_slow * 24.0 if math.isfinite(wall_days_slow) else float("nan"),
        wall_days_fast * 24.0 if math.isfinite(wall_days_fast) else float("nan"),
        atoms_estimate,
        box_edge,
        positive_count,
        negative_count,
    )


def build_plan(
    *,
    process_root: Path,
    microstate: str,
    sim_ns: float,
    anion_override: str | None = None,
    ion_count_override: int | None = None,
) -> tuple[list[RunEstimate], list[str]]:
    records = _load_records(process_root)
    if not records:
        raise ValueError(f"No historical simulations found under {process_root}")

    config = _load_microstate_config(process_root, microstate)
    if ion_count_override is not None:
        if ion_count_override <= 0:
            raise ValueError("--ion-count must be positive when provided.")
        config.positive_count = ion_count_override  # type: ignore[attr-defined]

    notes: list[str] = []
    results: list[RunEstimate] = []

    for ion_label, payload in ION_MATRIX.items():
        symbol = payload["symbol"]
        stats_records = _filter_by_ion(records, symbol)
        fallback_used = False
        if not stats_records:
            stats_records = records
            fallback_used = True
            notes.append(
                f"{ion_label}: no ion-specific runs found; using aggregate performance dataset."
            )

        density_stats = _density_stats(stats_records)
        perf_stats = _performance_stats(stats_records)

        cation_charge = _resolve_ion_charge(symbol, None, "cation")
        anion_symbol = (
            _normalize_ion_symbol(anion_override) if anion_override else config.anion
        )
        anion_charge = _resolve_ion_charge(anion_symbol, None, "anion")

        for conc_mM in payload["concentrations_mM"]:
            base_hours, slow_hours, fast_hours, atoms, box_edge, pos_count, neg_count = _compute_base_estimate(
                symbol,
                conc_mM,
                sim_ns=sim_ns,
                microstate_config=config,
                density_stats=density_stats,
                perf_stats=perf_stats,
                cation_charge=cation_charge,
                anion_symbol=anion_symbol,
                anion_charge=anion_charge,
            )
            for water_model, multiplier in WATER_MODELS:
                results.append(
                    RunEstimate(
                        ion_label=ion_label,
                        ion_symbol=symbol,
                        concentration_mM=conc_mM,
                        positive_count=pos_count,
                        negative_count=neg_count,
                        water_model=water_model,
                        multiplier=multiplier,
                        wall_hours=base_hours * multiplier,
                        wall_hours_slow=slow_hours * multiplier if math.isfinite(slow_hours) else float(
                            "nan"
                        ),
                        wall_hours_fast=fast_hours * multiplier if math.isfinite(fast_hours) else float(
                            "nan"
                        ),
                        atoms=atoms,
                        box_edge_nm=box_edge,
                        fallback_used=fallback_used,
                    )
                )

    return results, notes


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    default_process_root = repo_root / "ip6-run" / "process"

    parser = argparse.ArgumentParser(
        description=(
            "Assemble an 18-run simulation plan (3 ions × 3 concentrations × 2 water models) "
            "and estimate wall-clock times based on historical IP6 simulations."
        )
    )
    parser.add_argument(
        "--microstate",
        default="IP_010101",
        help="Reference microstate to source solute charge and default ion counts.",
    )
    parser.add_argument(
        "--sim-ns",
        type=float,
        default=100.0,
        help="Simulation length in nanoseconds for each run (default: 100).",
    )
    parser.add_argument(
        "--process-root",
        type=Path,
        default=default_process_root,
        help="Path to IP6 process directory with completed reference simulations.",
    )
    parser.add_argument(
        "--anion",
        help="Force a specific anion symbol (default: use microstate configuration).",
    )
    parser.add_argument(
        "--ion-count",
        type=int,
        help="Override positive ion count for box sizing (applies to all ions).",
    )
    args = parser.parse_args()

    results, notes = build_plan(
        process_root=args.process_root,
        microstate=args.microstate,
        sim_ns=args.sim_ns,
        anion_override=args.anion,
        ion_count_override=args.ion_count,
    )

    results.sort(key=lambda r: (r.ion_symbol, r.concentration_mM, r.water_model))

    header = (
        f"Simulation plan (microstate {args.microstate}, {args.sim_ns:.1f} ns each, T=310 K)"
    )
    print(header)
    print("-" * len(header))
    print(
        f"{'Ion':>6} {'Water':>12} {'Conc (mM)':>11} {'Cations':>9} {'Box edge (nm)':>14} "
        f"{'Atoms (~)':>12} {'Wall (typ)':>14} {'Wall (slow)':>14} {'Wall (fast)':>14}"
    )

    total_hours = 0.0
    total_hours_slow = 0.0
    total_hours_fast = 0.0
    total_hours_tip3p = 0.0
    total_hours_tip4p = 0.0

    for entry in results:
        total_hours += entry.wall_hours
        if math.isfinite(entry.wall_hours_slow):
            total_hours_slow += entry.wall_hours_slow
        if math.isfinite(entry.wall_hours_fast):
            total_hours_fast += entry.wall_hours_fast
        if entry.water_model == "TIP3P":
            total_hours_tip3p += entry.wall_hours
        else:
            total_hours_tip4p += entry.wall_hours

        atoms_fmt = _format_number(entry.atoms)
        slow_fmt = _format_hours(entry.wall_hours_slow) if math.isfinite(entry.wall_hours_slow) else "n/a"
        fast_fmt = _format_hours(entry.wall_hours_fast) if math.isfinite(entry.wall_hours_fast) else "n/a"
        fallback_marker = "*" if entry.fallback_used else " "

        print(
            f"{entry.ion_label:>6}{fallback_marker} "
            f"{entry.water_model:>12} "
            f"{entry.concentration_mM:11.3f} "
            f"{entry.positive_count:9d} "
            f"{entry.box_edge_nm:14.3f} "
            f"{atoms_fmt:>12} "
            f"{_format_hours(entry.wall_hours):>14} "
            f"{slow_fmt:>14} "
            f"{fast_fmt:>14}"
        )

    print()
    print(f"Total TIP3P wall-clock (median): {_format_hours(total_hours_tip3p)}")
    print(f"Total TIP4P/2005 wall-clock (median): {_format_hours(total_hours_tip4p)}")
    print(f"Combined total wall-clock (median): {_format_hours(total_hours)}")

    if total_hours_slow > 0:
        print(f"Combined slow-quartile estimate: {_format_hours(total_hours_slow)}")
    if total_hours_fast > 0:
        print(f"Combined fast-quartile estimate: {_format_hours(total_hours_fast)}")

    if notes:
        print("\nNotes:")
        for note in notes:
            print(f"- {note}")


if __name__ == "__main__":
    main()
