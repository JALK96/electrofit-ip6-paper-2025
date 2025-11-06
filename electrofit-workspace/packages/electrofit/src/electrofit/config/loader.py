
# packages/electrofit/src/electrofit/config/loader.py
from __future__ import annotations

import logging
from dataclasses import dataclass, field, is_dataclass, fields
from pathlib import Path
import typing as t

try:
    import tomllib  # py>=3.11
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore

"""TOML configuration loader (legacy .ef fallback removed)."""

if t.TYPE_CHECKING:  # pragma: no cover - typing only
    from electrofit.simulation.ions import DerivedSimulation

# -----------------
# Dataclass schema
# -----------------

@dataclass
class ProjectSection:
    name: str | None = None
    molecule_name: str | None = None
    residue_name: str | None = None
    charge: int | None = None
    protocol: str | None = None
    adjust_symmetry: bool = False
    ignore_symmetry: bool = False
    atom_type: str | None = None
    calculate_group_average: bool = False

@dataclass
class SymmetrySection:
    initial: str | None = None
    ensemble: str | None = None

@dataclass
class PathsSection:
    mdp_dir: str = "data/MDP"
    base_scratch_dir: str | None = None

# (HPC section removed for now)

# GROMACS runtime knobs (threads, pinning)
@dataclass
class GromacsRuntimeSection:
    threads: int = 16
    pin: bool = True

@dataclass
class GMXSection:
    runtime: GromacsRuntimeSection = field(default_factory=GromacsRuntimeSection)

# Simulation parameters (Step 3)
@dataclass
class SimulationBoxSection:
    type: str = "dodecahedron"
    edge_nm: float = 1.2

@dataclass
class IonTargetSection:
    desired_count: int | None = None
    concentration: float | None = None

@dataclass
class IonTargetsSection:
    positive_ion: IonTargetSection = field(default_factory=IonTargetSection)
    negative_ion: IonTargetSection = field(default_factory=IonTargetSection)

@dataclass
class SimulationIonsSection:
    cation: str = "NA"
    anion: str = "CL"
    salt_concentration: float | None = None
    enforce_neutrality: bool = True
    targets: IonTargetsSection | None = field(default_factory=IonTargetsSection)

@dataclass
class SimulationSection:
    box: SimulationBoxSection = field(default_factory=SimulationBoxSection)
    ions: SimulationIonsSection = field(default_factory=SimulationIonsSection)
    # Canonical location for forcefield selection (folder name with .ff)
    forcefield: str = "amber14sb.ff"
    derived: "DerivedSimulation | None" = None

@dataclass
class Config:
    project_root: Path
    project: ProjectSection = field(default_factory=ProjectSection)
    symmetry: SymmetrySection = field(default_factory=SymmetrySection)
    paths: PathsSection = field(default_factory=PathsSection)
    gmx: GMXSection = field(default_factory=GMXSection)
    simulation: SimulationSection = field(default_factory=SimulationSection)


# -----------------
# Helpers
# -----------------

def _load_toml(path: Path) -> dict:
    with path.open("rb") as f:
        return tomllib.load(f)


def _deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def _merge_into_dataclass(section, payload: dict):
    """Recursively merge a dict into a (possibly nested) dataclass instance."""
    for k, v in payload.items():
        if not hasattr(section, k):
            continue
        current = getattr(section, k)
        if is_dataclass(current) and isinstance(v, dict):
            _merge_into_dataclass(current, v)
        else:
            if v is not None:
                setattr(section, k, v)


def _flatten_dataclass(obj, prefix: str = ""):
    """Yield (key_path, value) for leaf attributes of nested dataclasses.

    Lists and non-dataclass objects are treated as leaves. Dataclass field order
    is preserved to keep dumps stable across runs (useful for diffing).
    """
    if is_dataclass(obj):
        for f in fields(obj):  # preserve declaration order
            val = getattr(obj, f.name)
            key = f"{prefix}.{f.name}" if prefix else f.name
            if is_dataclass(val):
                yield from _flatten_dataclass(val, key)
            else:
                yield key, val
    else:
        # Fallback – not expected for root call
        yield prefix or "value", obj


def dump_config(cfg: "Config", log_fn=print, header: bool = True):
    """Log all config settings (flattened) with a stable ordering.

    Format: [config] section.sub.key = value
    Nested dataclasses are traversed depth-first in field declaration order.
    """
    if header:
        log_fn("[config] -- begin full config dump --")
    for key, val in _flatten_dataclass(cfg):
        log_fn(f"[config] {key} = {val}")
    if header:
        log_fn("[config] -- end full config dump --")


def _infer_molecule_name_from_context(ctx: Path, project_root: Path) -> str | None:
    """
    Try to infer the molecule folder name from a context path located under either:
      - <project>/process/<mol>/...
      - <project>/data/input/<mol>/...
    Returns the molecule folder name if detected, else None.
    """
    try:
        ctx = ctx.resolve()
    except Exception:
        return None

    # Ensure we are inside the project tree
    try:
        ctx.relative_to(project_root)
    except Exception:
        # Not under project root; best-effort scan of parents
        pass

    # Walk up and look for 'process' or 'input' anchors
    p = ctx
    while True:
        if p.name == "process" and p.parent != p:
            # child of 'process' is the molecule name
            return ctx.relative_to(p).parts[0] if len(ctx.relative_to(p).parts) >= 1 else None
        if p.name == "input" and p.parent.name == "data":
            return ctx.relative_to(p).parts[0] if len(ctx.relative_to(p).parts) >= 1 else None
        if p.parent == p:
            break
        p = p.parent
    return None


# -----------------
# Loader
# -----------------

def load_config(
    project_root: t.Union[str, Path],
    config_path: t.Union[str, Path, None] = None,
    context_dir: t.Union[str, Path, None] = None,
    molecule_name: str | None = None,
) -> Config:
    root = Path(project_root).resolve()
    data: dict = {}
    compatibility_notes: list[str] = []

    ctx_path: Path | None = Path(context_dir).resolve() if context_dir else None
    mol_name: str | None = molecule_name
    if mol_name is None and ctx_path is not None:
        inferred = _infer_molecule_name_from_context(ctx_path, root)
        if inferred:
            mol_name = inferred

    project_toml = root / "electrofit.toml"
    provided = Path(config_path).resolve() if config_path else None

    tomls: list[Path] = []

    # 1) project-level defaults
    if project_toml.is_file():
        tomls.append(project_toml)

    # 2) data/input/<mol>/electrofit.toml
    if mol_name:
        input_mol_toml = root / "data" / "input" / mol_name / "electrofit.toml"
        if input_mol_toml.is_file():
            tomls.append(input_mol_toml)

    # 3) process/<mol>/electrofit.toml
    if mol_name:
        process_mol_toml = root / "process" / mol_name / "electrofit.toml"
        if process_mol_toml.is_file():
            tomls.append(process_mol_toml)

    # 4) per-run TOML in the context directory
    if ctx_path:
        run_toml = ctx_path / "electrofit.toml"
        if run_toml.is_file():
            tomls.append(run_toml)

    # 5) explicit --config (highest precedence)
    if provided and provided.is_file():
        tomls.append(provided)

    # De-duplicate while preserving order
    seen = set()
    ordered_unique: list[Path] = []
    for p in tomls:
        if p not in seen:
            ordered_unique.append(p)
            seen.add(p)

    logger = logging.getLogger(__name__)

    for p in ordered_unique:
        payload = _load_toml(p)
        # Compatibility shims are easier to handle before merging
        sim_payload = payload.get("simulation")
        if isinstance(sim_payload, dict):
            ions_payload = sim_payload.get("ions")
            if isinstance(ions_payload, dict):
                if "concentration" in ions_payload and "salt_concentration" not in ions_payload:
                    ions_payload["salt_concentration"] = ions_payload.pop("concentration")
                    compatibility_notes.append(
                        f"[config][compat] {p}: 'simulation.ions.concentration' is deprecated; use 'salt_concentration'."
                    )
        data = _deep_merge(data, payload)

    # If no TOML data found we simply return defaults (no legacy .ef fallback).

    # Build the Config object and merge sections (recursively for nested dataclasses)
    cfg = Config(project_root=root)

    # Accept both "gmx" and "gromacs" keys from TOML; merge "gromacs" into "gmx"
    if "gromacs" in data and "gmx" not in data:
        data["gmx"] = data.pop("gromacs")

    # Perform merges
    for section_name in ("project", "symmetry", "paths", "gmx", "simulation"):
        payload = data.get(section_name, {})
        section = getattr(cfg, section_name)
        if isinstance(payload, dict):
            _merge_into_dataclass(section, payload)
            setattr(cfg, section_name, section)

    if mol_name and not cfg.project.molecule_name:
        cfg.project.molecule_name = mol_name

    for note in compatibility_notes:
        logger.warning(note)

    # Derive simulation settings (ion counts, volumes, etc.)
    try:
        from electrofit.simulation.ions import derive_simulation_settings

        derived = derive_simulation_settings(
            cfg.simulation.box,
            cfg.simulation.ions,
            cfg.project.charge or 0,
        )
        cfg.simulation.derived = derived
        for line in derived.log_lines:
            logger.info(line)
        for warn in derived.warnings:
            logger.warning(warn)
    except ImportError:  # pragma: no cover - safety if module missing
        logger.debug("[config] ion derivation helper unavailable; skipping derived summary")
    except Exception as exc:
        raise ValueError(f"Ion configuration error: {exc}") from exc

    return cfg


__all__ = [
    "Config",
    "ProjectSection",
    "PathsSection",
    "SymmetrySection",
    "GMXSection",
    "SimulationSection",
    "load_config",
    "dump_config",
    "resolve_symmetry_flags",
]
_SYMMETRY_MODE_FLAGS: dict[str, tuple[bool, bool]] = {
    "antechamber": (False, False),
    "user": (True, False),
    "none": (True, True),
}


def _normalize_symmetry_mode(raw: str | None) -> str | None:
    if raw is None:
        return None
    mode = raw.strip().lower()
    return mode or None


def resolve_symmetry_flags(cfg: "Config", stage: t.Literal["initial", "ensemble"]) -> tuple[bool, bool]:
    """Resolve (adjust_symmetry, ignore_symmetry) for a pipeline stage.

    Falls back to legacy project-level booleans when the new [symmetry] section is absent.
    """
    section = getattr(cfg, "symmetry", None)
    mode_value = None
    if section is not None:
        mode_value = getattr(section, stage, None)
    mode = _normalize_symmetry_mode(mode_value)
    if mode:
        flags = _SYMMETRY_MODE_FLAGS.get(mode)
        if flags is None:
            raise ValueError(
                f"Invalid symmetry.{stage} mode '{mode_value}'. Expected one of: "
                + ", ".join(sorted(_SYMMETRY_MODE_FLAGS))
            )
        return flags
    # Legacy fallback: use project-level booleans
    return bool(getattr(cfg.project, "adjust_symmetry", False)), bool(getattr(cfg.project, "ignore_symmetry", False))
