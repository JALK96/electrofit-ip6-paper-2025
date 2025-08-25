# packages/electrofit/src/electrofit/config/loader.py
from __future__ import annotations

from dataclasses import dataclass, field, is_dataclass, fields
from pathlib import Path
import typing as t

try:
    import tomllib  # py>=3.11
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore

"""TOML configuration loader (legacy .ef fallback removed)."""

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
class SimulationIonsSection:
    cation: str = "NA"
    anion: str = "CL"
    concentration: float = 0.15

@dataclass
class SimulationSection:
    box: SimulationBoxSection = field(default_factory=SimulationBoxSection)
    ions: SimulationIonsSection = field(default_factory=SimulationIonsSection)
    # Canonical location for forcefield selection (folder name with .ff)
    forcefield: str = "amber14sb.ff"


@dataclass
class Config:
    project_root: Path
    project: ProjectSection = field(default_factory=ProjectSection)
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

    for p in ordered_unique:
        payload = _load_toml(p)
        data = _deep_merge(data, payload)

    # If no TOML data found we simply return defaults (no legacy .ef fallback).

    # Build the Config object and merge sections (recursively for nested dataclasses)
    cfg = Config(project_root=root)

    # Accept both "gmx" and "gromacs" keys from TOML; merge "gromacs" into "gmx"
    if "gromacs" in data and "gmx" not in data:
        data["gmx"] = data.pop("gromacs")

    # Perform merges
    for section_name in ("project", "paths", "gmx", "simulation"):
        payload = data.get(section_name, {})
        section = getattr(cfg, section_name)
        if isinstance(payload, dict):
            _merge_into_dataclass(section, payload)
            setattr(cfg, section_name, section)

    if mol_name and not cfg.project.molecule_name:
        cfg.project.molecule_name = mol_name

    return cfg