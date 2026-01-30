# Logging Refactor Plan (Steps 0–6)

This document is a design + execution plan for refactoring electrofit logging so that:

1) the per-step `process.log` files clearly answer *“what did my input do?”*,
2) “decisions” are not misleading (esp. symmetry during sampling),
3) logs are consistent across steps `step0`–`step6`,
4) we keep INFO logs concise (full config dumps remain DEBUG-only).

Scope: Step0–Step6 (Step7/8 exist, but are out of scope for this first refactor).

---

## 1) Problem statement (why this refactor)

Current pain points observed in real logs:

- **Step4 “symmetry” lines are misleading**: Step4 is sampling/extraction; it does not apply symmetry to charges, yet the decision summary prints something like `symmetry.effective: applied (antechamber defined)`.
- **Legacy booleans vs new config**: the project now uses
  - `[symmetry] initial = user|none|antechamber`
  - `[symmetry] ensemble = user|none|antechamber`
  but many logs still talk in terms of `project.adjust_symmetry` / `project.ignore_symmetry` booleans (which are derived flags and not what the user actually configured).
- **Users need a “story”**: per-step logs should make it obvious which config layers were used, which settings were effective, what artifacts were produced, and what is “planned for later”.

---

## 2) Current logging map (where logs live today)

### Project-level
- `<project>/step.log`
  - receives step summaries (and sometimes headers) from step scripts.

### Per molecule / per run directory

Step0:
- Writes into `<project>/step.log` only (no per-run folders yet).

Step1:
- Per molecule run dir: `process/<mol>/run_gau_create_gmx_in/process.log`
  - step1 orchestrator: `electrofit.pipeline.steps.step1`
  - domain work: `electrofit.domain.prep.process_initial` (includes DecisionModel logging)

Step2:
- Per molecule: `process/<mol>/run_gmx_simulation/process.log`
  - `electrofit.pipeline.steps.step2`

Step3:
- Per molecule: `process/<mol>/run_gmx_simulation/process.log`
  - `electrofit.pipeline.steps.step3`
  - GROMACS adapter: `electrofit.adapters.gromacs`

Step4:
- Per molecule: `process/<mol>/extracted_conforms/process.log`
  - `electrofit.pipeline.workers.step4_extract` (includes DecisionModel logging + cfg summary)
- Per conformer: directories under `process/<mol>/extracted_conforms/<conf>/` also get a snapshot, but Step4 mainly logs to the extracted_conforms root.

Step5:
- Per conformer: `process/<mol>/extracted_conforms/<conf>/process.log`
  - `electrofit.domain.charges.conformer_batch` sets up logging + snapshot
  - `electrofit.domain.charges.process_conformer` logs cfg + DecisionModel

Step6:
- Per molecule: `process/<mol>/results/process.log`
  - `electrofit.domain.aggregation.average_charges` logs cfg + DecisionModel

---

## 3) Proposed logging contract (applies to all steps)

### 3.1 Two concepts: `cfg` vs `decisions`

- `cfg` (configuration summary) = “effective values used by *this step*”
  - should prefer logging the *user-facing* knobs (e.g. symmetry modes `user|none|antechamber`) plus derived flags if helpful.
- `decisions` = “interpretation + consequences”
  - must distinguish:
    - **applied now** (this step actually uses it), vs
    - **planned later** (this step only prepares inputs for later steps).

### 3.2 Standard sections (INFO level)

Every per-step `process.log` should follow a consistent top-to-bottom structure:

1) **Header**
   - already exists: `electrofit <version> | step=<step> | git=<hash>`
2) **Context**
   - project root, run directory, molecule id, cwd, hostname (optional)
3) **Config provenance (snapshot sources)**
   - which files were used to seed + override `electrofit.toml` in this directory
4) **Effective cfg summary**
   - step-specific fields only (small, curated table)
5) **Decisions**
   - “what will happen / happened because of cfg”
6) **Actions & artifacts**
   - key actions taken, key files written/consumed (paths relative to run dir)
7) **End-of-step summary**
   - a single line suitable for grepping
8) **Warnings (as they occur)**
   - warnings should still be logged inline where they happen; optionally repeat a short “warning recap” at the end.

DEBUG level remains free to dump full config + tool commands.

---

## 4) Symmetry logging: required redesign

### 4.1 What users configure (preferred to log)

From `electrofit.toml`:
- `symmetry.initial = antechamber | user | none`
- `symmetry.ensemble = antechamber | user | none`

These should be logged directly, because they match user intent.

### 4.2 What the pipeline currently derives (still useful, but secondary)

The current resolver maps modes to legacy booleans:
- `antechamber -> (adjust=False, ignore=False)`
- `user -> (adjust=True, ignore=False)`
- `none -> (adjust=True, ignore=True)`

These derived flags can be logged *as derived output*, but they should not be the primary user-visible knob.

### 4.3 “Used here” vs “planned later”

We should log symmetry like:

- Step1: “symmetry.initial is **used** for initial charges”
- Step4: “symmetry.ensemble is **planned** (files copied / prepared), but not applied in sampling”
- Step5/6: “symmetry.ensemble is **used** for RESP / group averaging”

This eliminates the confusing “symmetry.effective: applied …” line in Step4.

### 4.4 DecisionModel changes (proposal)

Evolve `infra.decisions.DecisionModel` to add:

- `symmetry.initial.mode` (str | None)
- `symmetry.ensemble.mode` (str | None)
- `symmetry.stage_mode_used` (enum: `initial|ensemble|none`) or explicit booleans:
  - `symmetry.applied_in_this_step` (bool)
  - `symmetry.planned_for_later` (bool)
- optionally: `symmetry.source` (enum: `user_json | antechamber | missing`)

Then:
- Step4’s decision builder must set `symmetry.applied_in_this_step = False`
- Step4 can still emit `symmetry.planned_for_later = True` and warn if the required JSON is missing.

---

## 5) Step-by-step target logging content (Steps 0–6)

Below is what a *user* should be able to learn from each step’s logs, and which code should own it.

### Step0 (init workspace)

User questions to answer:
- Which molecules were discovered and copied?
- Was anything reused/skipped due to existing directories/locks?

Log content (INFO):
- Context: project root, `--molecule` filter if used
- Actions: copy `data/input/<mol> -> process/<mol>/run_gau_create_gmx_in`
- Outputs: list of created/updated dirs
- Summary: `copied=N reused=M skipped=K`

Code owners:
- `pipeline/steps/step0.py`

### Step1 (initial prep)

User questions:
- Which protocol ran (`bcc` vs `opt`) and what “charges origin” does that imply?
- Which symmetry mode was used for **initial** charges?
- Which external tools ran, and what artifacts were generated?

Log content (INFO):
- Config provenance: which electrofit.toml layers produced the snapshot
- Effective cfg: `protocol`, `net_charge`, `residue_name`, `atom_type`, `symmetry.initial.mode`
- Decisions: charges origin (`AM1-BCC` vs RESP), symmetry applied or suppressed
- Artifacts: `*.gcrt`, `*.gesp`, `*.esp`, `*.acpype/`, `*_GMX.{gro,itp,top}`, `posre_*.itp`, `initial_manifest.json`

Code owners:
- `pipeline/steps/step1.py` (snapshot + orchestrator)
- `domain/prep/process_initial.py` (DecisionModel + tool execution summary)

### Step2 (prepare gmx sim dir)

User questions:
- Which ACPYPE outputs were selected and copied?
- Which MDP templates were staged?

Log content (INFO):
- Effective cfg: MDP dir used, selected gro/top/itp/posre
- Artifacts: `run.json`, copied files list
- Summary per molecule: ready / skipped with reason

Code owners:
- `pipeline/steps/step2.py`

### Step3 (setup production MD)

User questions:
- Which forcefield and runtime settings were used?
- Were inputs found? If not, why skipped?
- What production directory/files were produced?

Log content (INFO):
- Effective cfg: forcefield, threads, pin, gpu flag, mdp_dir
- Actions: “set_up_production called”
- Summary: success/failure per molecule with reason

Code owners:
- `pipeline/steps/step3.py`
- `adapters/gromacs.py` (optional: emit a concise “what was generated” line)

### Step4 (extract conformers / sampling)

User questions:
- How were frames selected (linear/random/maxmin), how many, and with what seed?
- What is the total number of frames; which indices were picked?
- Which residue was extracted; did selection succeed?
- Which symmetry mode is **planned** for later RESP steps?

Log content (INFO):
- Effective cfg: sampling method/count/seed, residue_name (used for selection)
- Decisions:
  - sampling semantics (deterministic vs stochastic)
  - symmetry: `symmetry.ensemble.mode` planned (NOT applied here)
  - warn if needed files for later steps are missing (e.g. symmetry JSON absent but user asked for `user`)
- Actions: read `md_center.xtc` + `md.gro`, slice residue, sample indices
- Artifacts: list of conformer dirs created, plus each conformer’s snapshot presence

Code owners:
- `pipeline/workers/step4_extract.py`
- Decision builder: `infra/decisions.py` (`build_sampling_decision` must be corrected)

### Step5 (per-conformer Gaussian/RESP)

User questions:
- Which symmetry mode was actually applied during RESP?
- Did gaussian/espgen/resp succeed, and what was the runtime?
- What files should I inspect if something failed?

Log content (INFO):
- Effective cfg: conformer id, protocol, symmetry.ensemble.mode, net charge, residue name
- Decisions: symmetry applied or suppressed; charges origin
- Actions/artifacts: gaussian input/log, `.gesp`, `.esp`, resp outputs, updated mol2
- Fail-fast rule: if `.esp` is invalid, stop this conformer (no “0-charge” outputs)

Code owners:
- `domain/charges/process_conformer.py` (already logs cfg + DecisionModel)
- `cli/run_commands.py` for espgen validation/fail-fast behaviour (already implemented)

### Step6 (aggregate charges)

User questions:
- How many conformers contributed? Any missing/failed ones?
- What is the resulting charge set (avg, group-avg, cleaned)?
- Did group symmetry averaging happen (and why/why not)?

Log content (INFO):
- Effective cfg: protocol, symmetry.ensemble.mode, calculate_group_average, outlier removal flags
- Decisions: group averaging applied or skipped (missing JSON, no groups, etc.)
- Actions: read conformer outputs, compute means, optional outlier removal, optional group averaging
- Artifacts: `charges_dict.json`, `average_charges.chg`, `group_average_charges.chg`, `averaged_<mol>.mol2`, plots/histograms
- Summary: net charge before/after, number of conformers used

Code owners:
- `domain/aggregation/average_charges.py` (already logs cfg + DecisionModel; needs richer “inputs used / counts” summary)

---

## 6) Implementation roadmap (incremental, low-risk)

This is the suggested sequence to implement the redesign without breaking the pipeline.

### Phase A: define a stable log “shape”
1) Add a small helper in `infra/step_logging.py`:
   - `log_config_sources(step, snapshot_path, sources_used, overrides_used)`
   - (sources can be collected in `compose_snapshot` and returned as metadata)
2) Agree on curated `cfg` field lists per step (small tables).

### Phase B: symmetry + DecisionModel correctness
3) Extend `DecisionModel` to include symmetry mode strings (initial/ensemble) and applied/planned flags.
4) Fix `build_sampling_decision` to never claim symmetry is “applied”.
5) Update step4 extraction to log:
   - `symmetry.ensemble.mode` as “planned for step5/6”
   - warn when planned symmetry mode cannot be satisfied (e.g. mode=user but no JSON).

### Phase C: step-specific “action + artifacts” summaries
6) Step1: add a concise artifact summary line at the end of `process_initial`.
7) Step4: add frame count + selected indices summary (already partially present).
8) Step5: ensure conformer logs contain the key file names and “where to look if failed”.
9) Step6: add explicit counts (conformers used / missing) and net-charge normalization summary.

### Phase D: harmonize Step0–Step3
10) Standardize Step0–Step3 to emit:
    - config provenance,
    - curated cfg table,
    - end-of-step summary line.

---

## 7) Notes / non-goals

- Non-goal: dumping full config at INFO level (keep that DEBUG).
- Non-goal: changing scientific behaviour in this logging refactor (only semantics of logging and decision explanations).
- We should keep logs robust: logging failures must never crash the pipeline.

---

## 8) Keep existing “tool execution” logs (and improve later)

There are currently very useful low-level logs coming from the command wrappers / adapters, e.g.

- `New files/folders created: ...`
- `Executing command: respgen ...`
- `Executing python: electrofit.io.resp_edit.edit_resp_input(*args, **kwargs)`

**These logs should remain.** The Step0–Step6 logging refactor is primarily about:
- making `cfg` and `decisions` correct + readable,
- making provenance clear,
- and adding concise step summaries.

### Future improvements (de-prioritised)

Some “Executing python …” lines are not very informative because the arguments are not visible (or too opaque).
For a later iteration (not Phase A/B), we can enhance these logs by adding compact, user-readable context, e.g.:

- For file-to-file commands: log `input -> output` explicitly (filenames, not full paths).
- For `edit_resp_input`: log which RESP file was read, which was written, and a small “what changed” summary
  (e.g. changed N lines / toggled constraint block / added symmetry constraints), without spamming full diffs.
