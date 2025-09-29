from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Dict, Iterable

try:  # py311+
    import tomllib as _toml
except ModuleNotFoundError:  # pragma: no cover
    import tomli as _toml  # type: ignore

try:  # optional pretty TOML writer
    import tomli_w  # type: ignore
except ImportError:  # pragma: no cover
    tomli_w = None  # type: ignore

# Minimal TOML dump fallback (not full spec but enough for our nested dicts)
def _dump_toml(data: Dict[str, Any]) -> str:
    """Return TOML text for a nested dict.

    Prefer tomli_w if installed; otherwise use a lightweight deterministic
    writer that emits nested tables ([section] / [section.sub]). Only supports
    the scalar types we actually use (str, int, float, bool) and nested dicts.
    """
    if tomli_w:  # pragma: no cover - optional dependency path
        return tomli_w.dumps(data)  # type: ignore

    def _scalar(v: Any) -> str:
        if isinstance(v, bool):
            return "true" if v else "false"
        if isinstance(v, (int, float)):
            return str(v)
        # treat everything else as string
        s = str(v)
        s = s.replace("\\", "\\\\").replace("\"", "\\\"")
        return f'"{s}"'

    lines: list[str] = []

    # We emit tables depth-first; produce a list of (path_list, table_dict)
    def _walk(prefix: list[str], tbl: Dict[str, Any]):
        # Emit a header unless root (prefix empty) or this level has only nested dicts and no scalars? We'll still emit if prefix.
        scalars = [(k, v) for k, v in tbl.items() if not isinstance(v, dict)]
        if prefix:
            lines.append(f"[{'.'.join(prefix)}]")
        for k, v in sorted(scalars, key=lambda x: x[0]):
            lines.append(f"{k} = {_scalar(v)}")
        if scalars:
            lines.append("")
        # Recurse into dict children
        for k, v in sorted(tbl.items(), key=lambda x: x[0]):
            if isinstance(v, dict):
                _walk(prefix + [k], v)

    # Top-level: iterate keys in sorted order for determinism
    for key in sorted(data.keys()):
        val = data[key]
        if isinstance(val, dict):
            _walk([key], val)
        else:
            # Collect root scalars into a synthetic [root] section (unlikely in our schema)
            lines.append(f"{key} = {_scalar(val)}")

    # Remove potential trailing blank line
    while lines and lines[-1] == "":
        lines.pop()
    return "\n".join(lines) + "\n"


def _load(path: Path) -> Dict[str, Any]:
    with path.open('rb') as f:
        return _toml.load(f)  # type: ignore


def _walk(obj: Dict[str, Any], prefix: str = '') -> Iterable[tuple[str, Any]]:
    for k, v in obj.items():
        key = f'{prefix}.{k}' if prefix else k
        if isinstance(v, dict):
            yield from _walk(v, key)
        else:
            yield key, v


def _ensure_parent(dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    out = copy.deepcopy(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)  # type: ignore
        else:
            out[k] = copy.deepcopy(v)
    return out

RESTRICTED_DEFAULT = {
    'project.molecule_name',
    'project.charge',
    'simulation.ions.salt_concentration',
    'simulation.ions.targets.positive_ion.desired_count',
    'simulation.ions.targets.positive_ion.concentration',
    'simulation.ions.targets.negative_ion.desired_count',
    'simulation.ions.enforce_neutrality',
}


# INTERNAL: prefer using build_snapshot_with_layers; this direct merge kept for builder internals.
def merge_into_snapshot(
    snapshot_path: Path,
    override_path: Path | None,
    multi_molecule: bool,
    restricted_keys: set[str] | None = None,
    log_fn=print,
) -> bool:
    """Merge override TOML into existing snapshot file in-place.

    Returns True if any changes applied.
    """
    if not snapshot_path.is_file():
        return False
    if not override_path or not override_path.is_file():
        return False

    try:
        base = _load(snapshot_path)
    except Exception:
        base = {}
    try:
        over = _load(override_path)
    except Exception as e:  # pragma: no cover
        log_fn(f"[config][error] failed to read override {override_path}: {e}")
        return False

    restrict = restricted_keys or RESTRICTED_DEFAULT

    changes: list[tuple[str, Any, Any, bool]] = []  # key, old, new, restricted

    # Flatten override keys (no special aliasing; canonical forcefield lives under [simulation])
    def _apply(b: dict, o: dict, prefix: str = ''):
        for k, v in o.items():
            key = f'{prefix}.{k}' if prefix else k
            if isinstance(v, dict) and isinstance(b.get(k), dict):
                _apply(b[k], v, key)
            elif isinstance(v, dict) and k not in b:
                # whole subtree new
                b[k] = {}
                _apply(b[k], v, key)
            else:
                # scalar or list replace
                old = b.get(k)
                if old != v:
                    is_restricted = multi_molecule and key in restrict
                    changes.append((key, old, v, is_restricted))
                    b[k] = v

    # Work on a copy to detect net diff
    merged = _deep_merge(base, {})  # copy
    _apply(merged, over)

    if not changes:
        log_fn(f"[config] no overrides applied from {override_path}")
        return False

    # Persist
    _ensure_parent(snapshot_path)
    try:
        snapshot_path.write_text(_dump_toml(merged))
    except Exception as e:  # pragma: no cover
        log_fn(f"[config][error] failed writing merged snapshot: {e}")
        return False

    for key, old, new, restr in changes:
        if restr:
            log_fn(f"[config][warn] restricted key '{key}' override in multi-molecule context (old={old} -> new={new})")
        else:
            log_fn(f"[config][override] {key}: {old} -> {new}")
    return True


def fill_missing_into_snapshot(
    snapshot_path: Path,
    defaults_path: Path | None,
    log_fn=print,
) -> bool:
    """Fill ONLY missing keys in snapshot from defaults (no overrides).

    Returns True if any keys were added. Existing values are never changed.
    Nested dicts are created as needed. Logging uses [config][fill].
    """
    if not snapshot_path.is_file():
        return False
    if not defaults_path or not defaults_path.is_file():
        return False
    try:
        snap = _load(snapshot_path)
    except Exception:
        snap = {}
    try:
        defs = _load(defaults_path)
    except Exception as e:  # pragma: no cover
        log_fn(f"[config][error] failed to read defaults {defaults_path}: {e}")
        return False

    changes: list[tuple[str, Any]] = []

    def _fill(s: dict, d: dict, prefix: str = ''):
        for k, v in d.items():
            key = f'{prefix}.{k}' if prefix else k
            if isinstance(v, dict):
                if k not in s or not isinstance(s.get(k), dict):
                    s[k] = {}
                _fill(s[k], v, key)
            else:
                if k not in s:
                    s[k] = v
                    changes.append((key, v))

    _fill(snap, defs)
    if not changes:
        log_fn(f"[config] no fills applied from {defaults_path}")
        return False
    try:
        snapshot_path.write_text(_dump_toml(snap))
    except Exception as e:  # pragma: no cover
        log_fn(f"[config][error] failed writing filled snapshot: {e}")
        return False
    for key, val in changes:
        log_fn(f"[config][fill] {key}: {val}")
    return True
