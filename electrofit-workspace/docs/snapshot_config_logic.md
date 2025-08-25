# Snapshot Configuration Logic (Clean Reseed Behaviour)

Status: 2025-08-18 (updated after adoption of clean reseed semantics)
Scope: Documents the NEW deterministic snapshot layering (reseed every run) and preserves the legacy model for reference.

---
## New Core Mechanism: `compose_snapshot`
Location: `electrofit/infra/config_snapshot.py`

Algorithm (clean reseed):
1. Target snapshot path: `run_dir/electrofit.toml`.
2. ALWAYS determine a fresh seed each invocation from first existing source:
   `upstream -> molecule_input -> process_cfg -> project_defaults`.
   The previous local snapshot (if any) is ignored as seed material.
3. Apply STRONG overrides: `molecule_input -> process_cfg -> extra_override` (if present).
4. Apply FILL from `project_defaults` (adds missing keys/sections only).
5. Write snapshot (overwrite) and return path, or `None` if no seed source exists.

Properties:
- Upstream and defaults modifications propagate on every rerun automatically.
- Local manual edits to `run_dir/electrofit.toml` are ephemeral; persistence requires an override layer.
- Deterministic: rerunning without changing any seed/override/fill sources produces an identical file.

### Missing Source Handling
If `upstream` is specified but absent, the seed falls through to the next layer. If all layers are absent, snapshot creation fails (`None`).

### Override Precedence (lowest -> highest)
Seed < molecule_input < process_cfg < extra_override (CLI) ; fill adds only new keys from project defaults.

## Legacy Behaviour (for historical comparison)

Legacy algorithm (pre-change) used existing snapshot as seed if present, preventing upstream propagation after first creation. See final section for archived details.

Terminology:
- Base: The initial file contents after step 1–3.
- Strong override: Full value replacement merge (layer precedence above base).
- Fill: Adds missing keys only.

---
## Step-by-Step Behaviour (New Model)

For each step below:
- `RUN_DIR`: The directory whose snapshot is being (re)built.
- Files (possible sources):
  - `U` = upstream snapshot (previous step location as referenced by code)
  - `MI` = molecule input `data/input/<mol>/electrofit.toml`
  - `PC` = process level `process/<mol>/electrofit.toml`
  - `PD` = project defaults `<project_root>/electrofit.toml`
  - `EO` = `--config` extra override file (optional)
  - `RS` = existing `RUN_DIR/electrofit.toml`

Legend for scenario tables:
- `B=` base source chosen.
- `O:` applied strong overrides in order actually executed.
- `F:` fill phase.
- `Final:` resulting content conceptual order (earlier items may be shadowed by later overrides).

### Step0 (setup)
No snapshot composition occurs. Purely copies input molecule folders.

### Step1 (initial processing)
RUN_DIR: `process/<mol>/run_gau_create_gmx_in` (root per-molecule directory)
Parameters: `molecule_input=MI`, `project_defaults=PD`, optional `extra_override=EO`.
Seed selection: MI if present else PD. Each rerun repeats the same seed selection (ignores prior local snapshot content). Overrides & fill apply identically every run.

### Step2 (setup sim dir)
RUN_DIR: `process/<mol>/run_gmx_simulation`
Parameters: upstream (Step1 snapshot), MI, PC, PD, EO.
Each run seeds from current upstream if available; otherwise MI → PC → PD fallback. Upstream edits propagate automatically on rerun.

### Step3 (start sim)
Same RUN_DIR as Step2; clean reseed ensures any upstream change since Step2 is reflected immediately. The returned path may be logged; content is now meaningful (previously often redundant).

### Step4 (extract conformers)
RUN_DIR: `process/<mol>/extracted_conforms`
Parameters: upstream (simulation snapshot), MI, PC, PD, EO. Clean reseed ensures simulation changes flow into conformer extraction context immediately.

### Step5 (process conformers)
Per-conformer run directory: `process/<mol>/extracted_conforms/<conf>`.
`process_one` now composes a snapshot for each conformer directory before loading config (upstream = `extracted_conforms/electrofit.toml`). This harmonises layering with other steps; overrides and project defaults apply freshly per conformer on each invocation.

### Step6 (extract average charges)
RUN_DIR: `process/<mol>/results`
Parameters: `upstream=extracted_conforms/electrofit.toml`, `process_cfg=PC`, `molecule_input=MI`, `project_defaults=PD`, `extra_override=EO`.
Scenarios mirror Step2.

Special: Immediately after snapshot, `load_config(project_root, context_dir=results_dir, molecule_name=<mol>)` uses this snapshot as highest local context.

### Step7 (setup final sim)
RUN_DIR: `process/<mol>/run_final_gmx_simulation` (seed from results snapshot each run; overrides applied; defaults fill).

### Step8 (start final sim)
Clean reseed means upstream (results snapshot) changes propagate before final production simulation launch.

---
## Consolidated Comparison Table
| Step | Run Dir Snapshot Target | Upstream Source | Strong Overrides Order | Fill Source | Notes |
|------|-------------------------|-----------------|------------------------|-------------|-------|
| 0    | (n/a)                   | –               | –                      | –           | Pure copy only |
| 1    | run_gau_create_gmx_in   | –               | MI -> EO               | PD          | Root per-molecule seed |
| 2    | run_gmx_simulation      | Step1 snapshot  | MI -> PC -> EO         | PD          | Upstream seeds every run |
| 3    | run_gmx_simulation      | Step1 snapshot  | MI -> PC -> EO         | PD          | Reseed reflects upstream |
| 4    | extracted_conforms (inf)| Step2/3 snapshot | MI -> PC -> EO         | PD          | Inferred (not shown) |
| 5    | extracted_conformer dir | Step4 snapshot  | MI -> PC -> EO         | PD          | Per conformer (process_one) |
| 6    | results                 | Step4 snapshot  | MI -> PC -> EO         | PD          | Snapshot used for aggregation |
| 7    | run_final_gmx_simulation| Step6 results   | MI -> PC -> EO         | PD          | Prepares final inputs |
| 8    | run_final_gmx_simulation| Step6 results   | MI -> PC -> EO         | PD          | Launch final sim |

Legend: MI = molecule_input, PC = process_cfg, PD = project_defaults, EO = extra_override

---
## Remaining Considerations
1. Manual edits workflow changed (must migrate to override layers) – highlight in user docs.
2. Logging granularity for merges could be increased (future enhancement: structured diff logging).

---
## Illustrative Scenario Walkthroughs

### Scenario A: Add new key in project defaults after Step2, before rerunning Step3
Step3 reseeds from upstream then applies overrides and fill → new key appears (if absent upstream/overrides). Value mutations in project defaults also propagate because reseed no longer uses stale local snapshot.

### Scenario B: User adjusts a parameter in upstream (Step1 snapshot) after Step2 & Step3 ran
Re-run Step3: new upstream value present in seed → overrides applied → parameter updated (desired propagation).

### Scenario C: Providing `--config override.toml` in Step6
Overrides sequence ends with EO (strongest). Step7/8 reseed from results (upstream) and reapply EO if passed again. If EO not supplied later, its changes do not persist (expected – encourages explicit specification).

### Scenario D: Fresh molecule, full pipeline Step1→Step2→Step3
Each step reseeds from current upstream: Step2 from Step1 snapshot, Step3 again from Step1 snapshot (if unchanged) so typically identical unless upstream or overrides changed.

### Scenario E: Step5 processing with `--config` override
Each conformer directory reseeds snapshot; override file applied. Re-running with changed upstream (extracted_conforms snapshot) or defaults propagates automatically.

### Scenario F: Adding new key to `process/<mol>/electrofit.toml` after Step2, before Step3
- Step3 B=RS → override PC (adds/replaces key) → key updated (works because PC is in override sequence).

### Scenario G: Removing a key upstream intentionally
Rerun seeds from modified upstream lacking the key; if no override reintroduces it and defaults do not fill it, the key disappears (accurate reflection of upstream state).

---
## Summary of New Contract
* Snapshot is a transient, reproducible materialisation of declarative layers; it never serves as its own future seed.
* Upstream modifications propagate automatically.
* Override precedence: seed < MI < PC < CLI extra; PD fills only.
* Manual edits must move into MI / PC / CLI layer to persist.
* Step5 harmonised (per-conformer snapshots reseeded each run).

---
## Migration Guidance
1. Move any manual edits in run directory snapshots into `process/<mol>/electrofit.toml` or a CLI override file.
2. Re-run steps; verify new snapshot matches upstream expectations (diff if required).
3. Remove obsolete local snapshot commits from version control if they were previously tracked (should remain generated artifacts).

## Archived Legacy Description
The prior behaviour (sticky local snapshot, upstream ignored on reruns) is archived in the repository history (see pre-change version of this file) and summarised earlier. This section intentionally omits full repetition to reduce duplication.

---
*End of document.*
