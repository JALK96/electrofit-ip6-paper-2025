from __future__ import annotations

import argparse
import math
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

try:  # Python <3.11 support
    import tomllib  # type: ignore[attr-defined]
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore


AVOGADRO = 6.022_140_76e23  # mol^-1
SECONDS_PER_DAY = 86_400.0
CUBIC = "cubic"

ION_CHARGES = {
    "LI": 1,
    "NA": 1,
    "K": 1,
    "RB": 1,
    "CS": 1,
    "MG": 2,
    "CA": 2,
    "ZN": 2,
    "CL": -1,
    "BR": -1,
    "I": -1,
}


@dataclass
class SimulationRecord:
    microstate: str
    atoms: int
    volume_nm3: float
    cation: str
    positive_count: int
    concentration: float
    wall_time_s: float
    ns_per_day: float
    simulated_ns: float


@dataclass
class DensityStats:
    median: float
    q1: float
    q3: float


@dataclass
class PerformanceStats:
    median_k: float
    q1_k: float
    q3_k: float
    min_k: float
    max_k: float
    used_records: int
    target_ns: float


@dataclass
class MicrostateConfig:
    charge: int
    cation: str
    anion: str
    enforce_neutrality: bool
    positive_count: int
    negative_count: int | None
    box_type: str


def _normalize_ion_symbol(symbol: str) -> str:
    return "".join(ch for ch in symbol.upper() if ch.isalpha())


def _resolve_ion_charge(symbol: str, override: int | None, role: str) -> int:
    if override is not None:
        if override == 0:
            raise ValueError(f"{role.capitalize()} charge override must be non-zero.")
        return override
    if not symbol:
        raise ValueError(f"{role.capitalize()} symbol is required to infer a charge.")
    try:
        return ION_CHARGES[symbol]
    except KeyError as exc:
        raise ValueError(
            f"Unknown {role} '{symbol}'. Provide --{role.replace(' ', '-')}-charge explicitly."
        ) from exc


def _load_toml(path: Path) -> dict:
    with path.open("rb") as fh:
        return tomllib.load(fh)


def _quantile_triplet(values: Sequence[float]) -> tuple[float, float, float]:
    if not values:
        raise ValueError("Cannot compute quantiles of an empty collection.")
    if len(values) == 1:
        return values[0], values[0], values[0]
    q1, q2, q3 = statistics.quantiles(values, n=4, method="inclusive")
    return q1, q2, q3


def _find_tip3p_ions_gro(run_dir: Path) -> Path | None:
    for candidate in sorted(run_dir.glob("*_GMX_tip3p_ions.gro")):
        return candidate
    return None


def _parse_gro_summary(gro_path: Path) -> tuple[int, float]:
    with gro_path.open() as fh:
        lines = fh.readlines()
    if len(lines) < 3:
        raise ValueError(f"GRO file {gro_path} is too short.")
    atoms = int(lines[1].strip())
    try:
        box_vals = [float(x) for x in lines[-1].split()[:3]]
    except ValueError as exc:
        raise ValueError(f"Unable to parse box vectors in {gro_path}") from exc
    if len(box_vals) != 3:
        raise ValueError(f"Unexpected number of box vector entries in {gro_path}")
    volume = box_vals[0] * box_vals[1] * box_vals[2]
    return atoms, volume


def _parse_md_log(log_path: Path) -> tuple[float, float]:
    wall_time_s: float | None = None
    ns_per_day: float | None = None
    with log_path.open() as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("Time:"):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        wall_time_s = float(parts[2])
                    except ValueError:
                        continue
            elif stripped.startswith("Performance:"):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        ns_per_day = float(parts[1])
                    except ValueError:
                        continue
    if wall_time_s is None or ns_per_day is None:
        raise ValueError(f"Could not parse wall time or performance from {log_path}")
    return wall_time_s, ns_per_day


def _load_records(process_root: Path) -> list[SimulationRecord]:
    records: list[SimulationRecord] = []
    for micro_dir in sorted(process_root.glob("IP_*")):
        run_dir = micro_dir / "run_final_gmx_simulation"
        log_path = run_dir / "md.log"
        cfg_path = run_dir / "electrofit.toml"
        gro_path = _find_tip3p_ions_gro(run_dir)
        if not (log_path.exists() and cfg_path.exists() and gro_path and gro_path.exists()):
            continue
        try:
            cfg = _load_toml(cfg_path)
            atoms, volume = _parse_gro_summary(gro_path)
            wall_time_s, ns_per_day = _parse_md_log(log_path)
        except Exception:
            continue
        sim_section = cfg.get("simulation", {})
        ions_section = sim_section.get("ions", {})
        targets = ions_section.get("targets", {})
        positive_target = targets.get("positive_ion", {})
        record = SimulationRecord(
            microstate=micro_dir.name,
            atoms=atoms,
            volume_nm3=volume,
            cation=_normalize_ion_symbol(str(ions_section.get("cation", ""))),
            positive_count=int(positive_target.get("desired_count", 0) or 0),
            concentration=float(positive_target.get("concentration", 0.0) or 0.0),
            wall_time_s=wall_time_s,
            ns_per_day=ns_per_day,
            simulated_ns=ns_per_day * wall_time_s / SECONDS_PER_DAY,
        )
        records.append(record)
    return records


def _filter_by_ion(records: Iterable[SimulationRecord], ion: str) -> list[SimulationRecord]:
    normalized = _normalize_ion_symbol(ion)
    return [rec for rec in records if rec.cation == normalized]


def _density_stats(records: Iterable[SimulationRecord]) -> DensityStats:
    densities = [rec.atoms / rec.volume_nm3 for rec in records if rec.volume_nm3 > 0]
    if not densities:
        raise ValueError("No volume information available to compute density.")
    q1, _, q3 = _quantile_triplet(densities)
    return DensityStats(median=statistics.median(densities), q1=q1, q3=q3)


def _performance_stats(records: Iterable[SimulationRecord]) -> PerformanceStats:
    records = list(records)
    if not records:
        raise ValueError("No records available to compute performance statistics.")
    target_ns = statistics.median(rec.simulated_ns for rec in records)
    complete = [rec for rec in records if rec.simulated_ns >= 0.95 * target_ns and rec.ns_per_day > 0]
    if not complete:
        complete = [rec for rec in records if rec.ns_per_day > 0]
    k_values = [rec.ns_per_day * rec.atoms for rec in complete]
    if not k_values:
        raise ValueError("Performance dataset is empty.")
    q1, _, q3 = _quantile_triplet(k_values)
    return PerformanceStats(
        median_k=statistics.median(k_values),
        q1_k=q1,
        q3_k=q3,
        min_k=min(k_values),
        max_k=max(k_values),
        used_records=len(complete),
        target_ns=target_ns,
    )


def _load_microstate_config(process_root: Path, microstate: str) -> MicrostateConfig:
    cfg_path = process_root / microstate / "run_final_gmx_simulation" / "electrofit.toml"
    if not cfg_path.exists():
        raise FileNotFoundError(f"Missing electrofit.toml for microstate {microstate}")
    cfg = _load_toml(cfg_path)
    project_section = cfg.get("project", {})
    charge = int(project_section.get("charge", 0) or 0)
    sim_section = cfg.get("simulation", {})
    ions_section = sim_section.get("ions", {})
    box_section = sim_section.get("box", {})
    targets = ions_section.get("targets", {})
    positive_target = targets.get("positive_ion", {}) or {}
    negative_target = targets.get("negative_ion", {}) or {}

    positive_count = int(positive_target.get("desired_count", 0) or 0)
    if positive_count <= 0:
        raise ValueError(f"Positive ion count must be >0 for microstate {microstate}.")
    negative_count_raw = negative_target.get("desired_count")
    negative_count = int(negative_count_raw) if negative_count_raw is not None else None

    box_type = str(box_section.get("type", CUBIC) or CUBIC).lower()

    return MicrostateConfig(
        charge=charge,
        cation=_normalize_ion_symbol(str(ions_section.get("cation", ""))),
        anion=_normalize_ion_symbol(str(ions_section.get("anion", "CL")) or "CL"),
        enforce_neutrality=bool(ions_section.get("enforce_neutrality", True)),
        positive_count=positive_count,
        negative_count=negative_count,
        box_type=box_type,
    )


def _derive_box_lengths(volume_nm3: float, box_type: str) -> tuple[float, float, float]:
    if volume_nm3 <= 0:
        raise ValueError("Volume must be positive.")
    if box_type in {CUBIC, "cube"}:
        edge = volume_nm3 ** (1.0 / 3.0)
        return edge, edge, edge
    if box_type == "dodecahedron":
        # Volume of rhombic dodecahedron: (16/9)*sqrt(3)*a^3
        a = ((9.0 * volume_nm3) / (16.0 * math.sqrt(3.0))) ** (1.0 / 3.0)
        return a, a, a
    raise ValueError(f"Unsupported box type '{box_type}'.")


def _format_time(hours: float) -> str:
    if hours < 24:
        return f"{hours:.2f} h"
    days = hours / 24.0
    return f"{days:.2f} d"


def _print_dataset(records: Sequence[SimulationRecord]) -> None:
    header = f"{'microstate':>12} {'atoms':>8} {'vol[nm^3]':>10} {'ns/day':>10} {'wall[h]':>10}"
    print(header)
    for rec in records:
        wall_hours = rec.wall_time_s / 3600.0
        print(
            f"{rec.microstate:>12} {rec.atoms:8d} {rec.volume_nm3:10.3f} "
            f"{rec.ns_per_day:10.3f} {wall_hours:10.2f}"
        )


def estimate_runtime(
    ion_type: str,
    concentration: float,
    microstate: str,
    sim_ns: float,
    process_root: Path,
    override_positive: int | None = None,
    show_dataset: bool = False,
    *,
    anion_type: str | None = None,
    cation_charge: int | None = None,
    anion_charge: int | None = None,
) -> None:
    if concentration <= 0:
        raise ValueError("Concentration must be positive (mol/L).")
    if sim_ns <= 0:
        raise ValueError("Simulation length must be positive (ns).")

    records = _load_records(process_root)
    if not records:
        raise ValueError(f"No historical records found under {process_root}.")

    normalized_cation = _normalize_ion_symbol(ion_type)
    if not normalized_cation:
        raise ValueError("Ion type must contain alphabetic characters (e.g. 'Na', 'Ca').")

    ion_specific_records = _filter_by_ion(records, normalized_cation)
    fallback_used = False
    stats_records = ion_specific_records
    if not stats_records:
        stats_records = records
        fallback_used = True

    if show_dataset:
        label = "Historical dataset used for scaling"
        if fallback_used:
            label += " (fallback to all available cations)"
        print(f"{label}:")
        _print_dataset(stats_records)
        print()

    density_stats = _density_stats(stats_records)
    perf_stats = _performance_stats(stats_records)

    config = _load_microstate_config(process_root, microstate)

    cation_charge_val = _resolve_ion_charge(normalized_cation, cation_charge, "cation")
    anion_symbol = _normalize_ion_symbol(anion_type) if anion_type else config.anion
    if not anion_symbol:
        raise ValueError("Anion symbol could not be determined; specify --anion.")
    anion_charge_val = _resolve_ion_charge(anion_symbol, anion_charge, "anion")
    if anion_charge_val >= 0:
        raise ValueError("Anion charge must be negative (e.g. -1 for Cl-).")
    if cation_charge_val <= 0:
        raise ValueError("Cation charge must be positive (e.g. +1 for Na+, +2 for Ca2+).")

    positive_count = override_positive if override_positive is not None else config.positive_count
    adjustments: list[str] = []

    if config.enforce_neutrality:
        min_positive = math.ceil(-config.charge / cation_charge_val) if cation_charge_val > 0 else 0
        if positive_count < min_positive:
            adjustments.append(
                f"Increased positive ion count to {min_positive} to satisfy neutrality "
                f"with solute charge {config.charge:+d} e."
            )
            positive_count = min_positive

        modulus = abs(anion_charge_val)
        neutrality_term = positive_count * cation_charge_val + config.charge
        remainder = neutrality_term % modulus
        if remainder != 0:
            for delta in range(1, modulus + 1):
                test_positive = positive_count + delta
                test_term = test_positive * cation_charge_val + config.charge
                if test_term % modulus == 0:
                    positive_count = test_positive
                    neutrality_term = test_term
                    adjustments.append(
                        f"Adjusted positive ion count by +{delta} to make anion count integral."
                    )
                    break
            else:
                raise ValueError(
                    "Unable to satisfy neutrality with integer ion counts. "
                    "Override --ion-count or adjust ion charges."
                )
        total_positive_charge = positive_count * cation_charge_val
        if total_positive_charge + config.charge < 0:
            raise ValueError(
                "Positive ion count insufficient to offset solute charge; "
                "increase --ion-count or choose a higher concentration."
            )
        negative_count = (total_positive_charge + config.charge) // abs(anion_charge_val)
        net_charge = total_positive_charge + config.charge + negative_count * anion_charge_val
    else:
        negative_count = config.negative_count or 0
        total_positive_charge = positive_count * cation_charge_val
        net_charge = total_positive_charge + config.charge + negative_count * anion_charge_val

    mol_positive = positive_count / AVOGADRO
    volume_litre = mol_positive / concentration
    volume_nm3 = volume_litre * 1e24
    lengths_nm = _derive_box_lengths(volume_nm3, config.box_type)

    atoms_median = density_stats.median * volume_nm3
    atoms_q1 = density_stats.q1 * volume_nm3
    atoms_q3 = density_stats.q3 * volume_nm3

    ns_day_median = perf_stats.median_k / atoms_median
    ns_day_slow = perf_stats.q1_k / atoms_q3 if atoms_q3 > 0 else float("nan")
    ns_day_fast = perf_stats.q3_k / atoms_q1 if atoms_q1 > 0 else float("nan")
    ns_day_worst = perf_stats.min_k / atoms_q3 if atoms_q3 > 0 else float("nan")
    ns_day_best = perf_stats.max_k / atoms_q1 if atoms_q1 > 0 else float("nan")

    wall_days_median = sim_ns / ns_day_median
    wall_days_slow = sim_ns / ns_day_slow if ns_day_slow > 0 else float("nan")
    wall_days_fast = sim_ns / ns_day_fast if ns_day_fast > 0 else float("nan")
    wall_days_worst = sim_ns / ns_day_worst if ns_day_worst > 0 else float("nan")
    wall_days_best = sim_ns / ns_day_best if ns_day_best > 0 else float("nan")

    print(f"Microstate: {microstate}")
    print(f"  Solute charge: {config.charge:+d} e")
    print(f"  Ion type: {normalized_cation} (historical records: {len(stats_records)})")
    if fallback_used:
        print("  Note: No ion-specific runs found; using aggregate dataset across microstates.")
    print(f"  Cation charge: {cation_charge_val:+d} e")
    print(f"  Anion: {anion_symbol} ({anion_charge_val:+d} e)")
    print(f"  Neutrality enforced: {config.enforce_neutrality}")
    if adjustments:
        for note in adjustments:
            print(f"  Adjustment: {note}")
    print(
        f"  Positive ions: {positive_count} "
        f"(override={'yes' if override_positive is not None else 'no'})"
    )
    print(f"  Negative ions: {negative_count}")
    print(f"  Net system charge estimate: {net_charge:+d} e")
    print(f"  Target concentration: {concentration:.3f} mol/L")

    print()
    print("Derived box (from ion targets):")
    print(f"  Volume: {volume_nm3:.3f} nm^3 ({volume_litre:.3e} L)")
    print(f"  Box type: {config.box_type}")
    print(f"  Box lengths (nm): {', '.join(f'{x:.3f}' for x in lengths_nm)}")

    print()
    print("Estimated composition:")
    print(f"  Atoms (median density): {int(round(atoms_median))}")
    print(
        f"  Atoms IQR: {int(round(atoms_q1))} – {int(round(atoms_q3))} "
        "(scaled from historical densities)"
    )

    print()
    print(f"Predicted runtime for {sim_ns:.1f} ns (based on historical wall-clock performance):")
    print(
        f"  Wall clock (typical): {_format_time(wall_days_median * 24)} "
        f"(~{ns_day_median:.1f} ns/day)"
    )
    if math.isfinite(wall_days_slow):
        print(
            f"  Wall clock (slow quartile): {_format_time(wall_days_slow * 24)} "
            f"(~{ns_day_slow:.1f} ns/day)"
        )
    if math.isfinite(wall_days_fast):
        print(
            f"  Wall clock (fast quartile): {_format_time(wall_days_fast * 24)} "
            f"(~{ns_day_fast:.1f} ns/day)"
        )
    if math.isfinite(wall_days_worst):
        print(
            f"  Wall clock (worst observed): {_format_time(wall_days_worst * 24)} "
            f"(~{ns_day_worst:.1f} ns/day)"
        )
    if math.isfinite(wall_days_best):
        print(
            f"  Wall clock (best observed): {_format_time(wall_days_best * 24)} "
            f"(~{ns_day_best:.1f} ns/day)"
        )

    dataset_label = (
        f"{len(stats_records)} run(s) across microstates"
        if fallback_used
        else f"{len(stats_records)} {normalized_cation}-specific run(s)"
    )
    print()
    print(
        f"Historical dataset used: {perf_stats.used_records} complete runs "
        f"(source: {dataset_label}; median simulated length {perf_stats.target_ns:.1f} ns)."
    )


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    default_process_root = repo_root / "ip6-run" / "process"

    parser = argparse.ArgumentParser(
        description=(
            "Estimate simulation wall-clock time for IP6 microstates by "
            "scaling historical runs with different ion concentrations."
        )
    )
    parser.add_argument("ion_type", help="Cation symbol (e.g. Na, Ca, Mg).")
    parser.add_argument("concentration", type=float, help="Desired molarity (mol/L).")
    parser.add_argument(
        "--microstate",
        default="IP_010101",
        help="Reference microstate to pull solute charge and default ion count from.",
    )
    parser.add_argument(
        "--sim-ns",
        type=float,
        default=100.0,
        help="Simulation length in nanoseconds to estimate (default: 100).",
    )
    parser.add_argument(
        "--ion-count",
        type=int,
        help="Override positive ion count instead of using the microstate default.",
    )
    parser.add_argument(
        "--anion",
        help="Override anion symbol (default: taken from the microstate configuration).",
    )
    parser.add_argument(
        "--cation-charge",
        type=int,
        help="Override cation charge in elementary charge units (e).",
    )
    parser.add_argument(
        "--anion-charge",
        type=int,
        help="Override anion charge in elementary charge units (e).",
    )
    parser.add_argument(
        "--process-root",
        type=Path,
        default=default_process_root,
        help="Path to the IP6 process directory (default: repo/ip6-run/process).",
    )
    parser.add_argument(
        "--show-dataset",
        action="store_true",
        help="Print the historical dataset used for the estimate.",
    )
    args = parser.parse_args()

    estimate_runtime(
        ion_type=args.ion_type,
        concentration=args.concentration,
        microstate=args.microstate,
        sim_ns=args.sim_ns,
        process_root=args.process_root,
        override_positive=args.ion_count,
        show_dataset=args.show_dataset,
        anion_type=args.anion,
        cation_charge=args.cation_charge,
        anion_charge=args.anion_charge,
    )


if __name__ == "__main__":
    main()
