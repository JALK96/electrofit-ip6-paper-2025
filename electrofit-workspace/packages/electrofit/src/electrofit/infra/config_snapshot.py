"""Config Snapshot Composition (infrastructure layer).

Public entry:
    compose_snapshot(...) -> Path|None
    CONFIG_ARG_HELP (Reuse in CLI help strings)

2025-08 (behaviour change): Reruns now *always* re-seed the snapshot from the
highest-precedence available *source* (upstream > molecule_input > process_cfg >
project_defaults) instead of preserving an existing local snapshot file as
immutable base. This yields a "clean" deterministic override each run. The
previous local snapshot is ignored (overwritten) unless **no** source layer
exists (fallback).
"""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Callable
from electrofit.config.merger import merge_into_snapshot, fill_missing_into_snapshot

CONFIG_ARG_HELP = (
    "Pfad zu einer zusätzlichen electrofit.toml: Molekül- und Prozess-Dateien überschreiben Projekt-Defaults; "
    "diese Datei wirkt als zusätzliche stärkste Override-Ebene (ändert bestehende Werte, füllt nichts nur auf)."
)

# Internal helper

def _seed_snapshot(target_dir: Path, sources: Iterable[Path]) -> Path | None:
    """(Re)create snapshot from first existing source.

    Overwrites existing snapshot unconditionally (clean behaviour) except when
    *no* source exists; then keeps existing snapshot if present, else returns None.
    """
    snap = target_dir / "electrofit.toml"
    for cand in sources:
        if cand and cand.is_file():
            try:
                snap.write_text(cand.read_text())
            except Exception as e:  # log minimal diagnostic via print (no logger guaranteed here)
                print(f"[snapshot][warn] failed to seed from {cand}: {e}")
            break
    else:  # no source found
        if not snap.exists():
            return None
    return snap


def compose_snapshot(
    run_dir: Path,
    project_root: Path,
    molecule: str | None,
    multi_molecule: bool,
    log_fn: Callable[[str], None],
    upstream: Path | None = None,
    process_cfg: Path | None = None,
    molecule_input: Path | None = None,
    project_defaults: Path | None = None,
    extra_override: Path | None = None,
) -> Path | None:
    """Construct (or re-construct) snapshot in ``run_dir``.

    Clean semantics (idempotent but *not* sticky): Every invocation reseeds from
    first existing source among: upstream, molecule_input, process_cfg, project_defaults.
    Existing local snapshot content is ignored (overwritten) so upstream & lower
    layers propagate on reruns.

    Strong override sequence:
        molecule_input -> process_cfg -> extra_override
    Fill phase (adds only missing keys):
        project_defaults
    """
    sources: list[Path] = []
    if upstream:
        sources.append(upstream)
    if molecule_input:
        sources.append(molecule_input)
    if process_cfg:
        sources.append(process_cfg)
    if project_defaults:
        sources.append(project_defaults)
    snap = _seed_snapshot(run_dir, sources)
    if not snap or not snap.is_file():
        return None
    # Strong overrides
    for layer in [molecule_input, process_cfg, extra_override]:
        if layer and layer.is_file():
            try:
                merge_into_snapshot(snap, layer, multi_molecule=multi_molecule, log_fn=log_fn)
            except Exception as e:
                log_fn(f"[snapshot][warn] merge failed for {layer}: {e}")
    # Fill missing
    if project_defaults and project_defaults.is_file():
        try:
            fill_missing_into_snapshot(snap, project_defaults, log_fn=log_fn)
        except Exception as e:
            log_fn(f"[snapshot][warn] fill failed: {e}")
    return snap

__all__ = ["compose_snapshot", "CONFIG_ARG_HELP"]
