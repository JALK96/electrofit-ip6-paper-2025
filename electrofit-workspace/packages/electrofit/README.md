# electrofit (core)

Minimal tooling to parameterise small molecules/residues, run short GROMACS MD, extract conformers, and drive RESP charge aggregation.

## What lives in this package
- `cli/`: the `electrofit` command dispatcher (step0–step8, REMD variants) and helper CLIs for symmetry / charge editing and plotting.
- `config/`: loading + layering of `electrofit.toml`, derived ion/box settings, snapshot helpers.
- `pipeline/`: step orchestration (prep, MD, conformer extraction, RESP, aggregation, final MD, REMD).
- `domain/`: pure logic (sampling, symmetry constraints, charge processing, final/remd sim manifests).
- `adapters/`: thin wrappers over external tools (GROMACS, Gaussian, RESP) and result helpers.
- `io/`: MOL2/RESP/equiv group readers & writers, topology utilities.
- `infra/`: logging, scratch management, decisions/snapshots.
- `plotting/` and `viz/`: optional plotting utilities for charge histograms and helper glyphs.

## Installation (editable example)
```bash
pip install -e ./packages/electrofit
```

## CLI workflow (high level)
Every command requires `--project <project_root>` and accepts an optional `--config <path>` override. Molecule name is inferred from context (`process/<mol>` or `data/input/<mol>`) when possible.

| Command      | Purpose (inputs ➜ outputs) |
|--------------|----------------------------|
| `step0`      | Scaffold `process/` tree for each molecule under `data/input/`; idempotent copy/normalise. |
| `step1`      | Initial prep (Gaussian/RESP staging, symmetry JSON, ACPYPE-ready files). Outputs `process/<mol>/run_gau_create_gmx_in/`. |
| `step2`      | Build `run_gmx_simulation` directories: stage MDPs from `[paths].mdp_dir`, create run manifest, copy topology/start structures. |
| `step3`      | EM → NVT → NPT → short production MD. Respects `[gmx.runtime]` threads/pin/gpu. Produces `md.gro`, `md_center.xtc`. |
| `step4`      | Sample conformers from production trajectory into `extracted_conforms/` using `linear|random|maxmin` frame selection. |
| `step5`      | Per-conformer Gaussian/RESP pipeline (parallelisable), writes RESP charges per conformer. |
| `step6`      | Aggregate charges (mean, optional symmetry/group averaging, optional outlier removal), update MOL2/ACPYPE. |
| `step7`      | Final GROMACS setup using aggregated charges (`run_final_gmx_simulation/`). |
| `step8`      | Final MD production with aggregated charges. |
| `step7remd`  | REMD setup (templates, replica folders, `remd*.tpr`). |
| `step8remd`  | REMD equilibration + prep, writes `run_remd.sh`; can launch via generated script. |

### Minimal workflow example
```bash
# 1) scaffold and prep
electrofit step0 --project /path/to/project
electrofit step1 --project /path/to/project

# 2) build and run short MD, then sample conformers
electrofit step2 --project /path/to/project
electrofit step3 --project /path/to/project
electrofit step4 --project /path/to/project --sampling-method maxmin --sample 25 --seed 123

# 3) RESP ensemble + aggregation
electrofit step5 --project /path/to/project --workers 8
electrofit step6 --project /path/to/project --plot-histograms --remove-outlier
```

## Configuration reference (`electrofit.toml`)
Configuration is layered (lowest → highest precedence):
1. `<project>/electrofit.toml`
2. `<project>/data/input/<mol>/electrofit.toml`
3. `<project>/process/<mol>/electrofit.toml`
4. `<run_dir>/electrofit.toml` (if present)
5. `--config <path>`

Fields are merged into nested dataclasses; unknown keys are ignored. Defaults shown in parentheses.

### [project]
- `name` (None): free text project label.
- `molecule_name` (None): overrides inferred molecule folder name.
- `residue_name` (None): residue name used in topologies.
- `charge` (None): total solute charge (integer).
- `protocol` (None): charge protocol tag (`bcc`, `opt`, etc.).
- `atom_type` (None): optional forcefield atom type hint.
- `adjust_symmetry` (false), `ignore_symmetry` (false): legacy symmetry flags (overridden by `[symmetry]`).
- `calculate_group_average` (false): enable symmetry-group averaging during aggregation.

### [symmetry]
- `initial` (`antechamber`|`user`|`none`, default None → falls back to project flags): mode for the initial (step1) stage.
- `ensemble` (`antechamber`|`user`|`none`, default None): mode for ensemble/aggregation stage. Modes map to `(adjust_symmetry, ignore_symmetry)` as follows: `antechamber` → (False, False); `user` → (True, False); `none` → (True, True).

### [paths]
- `mdp_dir` (`"data/MDP"`): directory containing MDP templates relative to project root.
- `base_scratch_dir` (None): root for scratch working dirs (env vars expanded).

### [gmx.runtime]
- `threads` (16): OpenMP threads per `gmx_mpi` rank.
- `pin` (true): CPU pinning flag passed to GROMACS.
- `gpu` (false): whether to request GPU aware settings in templates (templates must support it).

### [simulation.box]
- `type` (`"dodecahedron"`): supports `dodecahedron` or `cubic` for automatic sizing.
- `edge_nm` (1.2): padding distance when computing box size around solute.

### [simulation.ions]
- `cation` (`"NA"`), `anion` (`"CL"`).
- `salt_concentration` (default 0.15 M if unset): bulk salt mode (legacy default when no targets are provided).
- `enforce_neutrality` (true).
- `targets`: enables target-count mode if any target has payload.
  - `positive_ion.desired_count` (int, required in target mode)
  - `positive_ion.concentration` (float, mol/L, required in target mode)
  - `negative_ion.desired_count` (int, optional; inferred if neutrality enforced)
  - `negative_ion.concentration` (not supported; leave unset)
- Derived summary (ion counts, inferred neutrality, effective molarity, box lengths) is logged via `derive_simulation_settings`.

### [simulation]
- `forcefield` (`"amber14sb.ff"`): folder name under `share/gromacs/top` or project-supplied.

### [remd]
- `enabled` (false)
- `nreplicas` (16)
- `tmin` (300.0), `tmax` (420.0)
- `spacing` (`"geometric"` or `"list"`)
- `temperatures` (list, used when `spacing="list"`)
- `replex` (100): exchange period (steps)
- `mdp_template` (`"REMD.mdp"`)
- `mpi_ranks` (None): overrides `nreplicas` when launching mpirun for REMD.

### [sampling] (used by step4 only)
- `method` (`"linear"` | `"random"` | `"maxmin"`, default `linear`)
- `count` (int, default 20)
- `seed` (int | None): seed for random/maxmin; `None` → deterministic linear or library default.

## Sampling quick reference (step4)
```bash
electrofit step4 --project <proj> \
  --sampling-method maxmin \
  --sample 30 \
  --seed 42 \
  --workers 4
```
Outputs PDBs under `process/<mol>/extracted_conforms/` and logs a sampling decision record (method, count, seed, symmetry notes).

## Notes on configuration resolution
- Unknown keys are ignored (safe to keep project-level extras).
- CLI flags always win over TOML.
- Snapshot composition for per-run folders follows: upstream → molecule input → process cfg → project defaults → explicit override; existing snapshots are re-seeded on each run (clean behaviour).

## Minimal API surface (Python)
```python
from electrofit.config.loader import load_config, dump_config

cfg = load_config(project_root="/path/to/project")
dump_config(cfg)  # prints flattened keys and values

from electrofit.simulation.ions import derive_simulation_settings
derived = derive_simulation_settings(cfg.simulation.box, cfg.simulation.ions, cfg.project.charge or 0)
print(derived.ions.positive_count, derived.ions.negative_count)
```

The public entrypoint `python -m electrofit` is equivalent to running the `electrofit` console script.

## Logging and scratch
- Default logging is INFO-level stdout unless a handler is already configured by the caller.
- Scratch directories are managed by `infra.scratch_manager` under `[paths].base_scratch_dir` (env-expanded); operations are idempotent.

## Safety and compatibility
- Legacy `.ef` configs are no longer read; use TOML only.
- Deprecated key shim: `simulation.ions.concentration` is auto-mapped to `simulation.ions.salt_concentration` with a warning.
- Symmetry defaults fall back to project-level booleans when `[symmetry]` is omitted.

