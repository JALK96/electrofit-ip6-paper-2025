# ip6-analysis Refactor Plan (Correctness-First)

## 1) Scope and Goal

This plan focuses on:

1. Fixing the H-bond XPM/log row-alignment bug everywhere it affects analysis.
2. Adding a `terminal-only` mode for P->P matrices.
3. Removing duplicated H-bond parsing logic to prevent future divergence.
4. Refactoring in safe, testable increments.

Primary objective: **all generated P->P and row-indexed H-bond outputs must be physically consistent and reproducible across microstates**.

---

## 2) Double-Check Findings (Plan Risks and Corrections)

### Confirmed risks in the original plan

1. **Global flip risk**  
   Blindly reversing XPM rows everywhere could break workflows if a parser or upstream source already compensates orientation.

2. **Hidden dependency risk**  
   Multiple modules duplicate `parse_xpm`, `parse_hbond_log_to_dataframe`, `analyze_hydrogen_bonds`. Fixing only one path leaves silent inconsistencies.

3. **Mapping semantics risk**  
   `terminal-only` selection can drift if oxygen classification is hardcoded inconsistently across modules or molecule variants.

4. **Regression risk**  
   No test suite currently exists in package; large refactor without tests is high risk.

5. **Operational risk**  
   Environment-level `gmx` runtime is currently not guaranteed (`libmpi_cxx.so.20` issue seen), so validation must not depend on rerunning full GROMACS pipelines.

### Plan corrections after double-check

1. Do not keep any “old/wrong row assignment” runtime mode.
2. Introduce a **single canonical H-bond IO utility** and migrate all callers.
3. Add fixture-based tests for row alignment + mapping before broad cleanup.
4. Separate “correctness changes” from “style/cleanup changes” into phases.
5. Validate using existing produced files (log/xpm/idx) and deterministic recomputation where possible.

---

## 3) Design Decisions

1. **Canonical parser module**
   - New shared module under `src/electrofit_analysis/structure/util/` (e.g. `hbond_io.py`).
   - Contains:
     - robust XPM parsing
     - canonical row orientation handling
     - hb.log parsing utilities
     - optional matrix helper utilities for row-indexed workflows

2. **Row orientation policy**
   - Canonical parser returns rows aligned to `hb.log` / `hb_idx` order.
   - No user-facing switch for legacy orientation.
   - Add validation helper(s) in tests to prevent reintroduction.

3. **Terminal-only mode**
   - Add explicit oxygen selection mode for P->P matrix:
     - `all` (default): bridging + terminal oxygens
     - `terminal`: only terminal oxygens
   - Implement in `pp-matrix` command path only for now.

4. **Incremental migration**
   - Keep temporary thin wrappers in old modules while migrating imports.
   - Remove wrappers only after tests and callers are fully updated.

---

## 4) Phase Plan

## Phase A - Guardrails (tests + fixtures) [BLOCKER BEFORE MASS EDITS]

- [x] A1. Create `tests/` package and baseline test config.
- [x] A2. Add compact fixture files for H-bond parser tests (synthetic + one real microstate subset if lightweight).
- [x] A3. Add tests:
  - [x] XPM parser shape/color handling.
  - [x] Row-order alignment against log order.
  - [ ] `hbonds_per_index` row correspondence.
  - [x] P->P mapping smoke test for `IP_101101` behavior (specifically non-high `P4->P2` after fix).
- [x] A4. Add tests for oxygen mode:
  - [x] `all` vs `terminal` yields expected differences.

Acceptance criteria for Phase A:

- At least one failing test demonstrates current bug before fix.
- Tests can run without GROMACS binaries.

---

## Phase B - Correctness Implementation (H-bond stack)

- [x] B1. Add canonical shared module `hbond_io.py`.
- [x] B2. Implement canonical XPM parser with correct row alignment semantics.
- [x] B3. Move/centralize shared parsing logic now duplicated in:
  - `cli/h_bonds/h_bonds_ip6.py`
  - `structure/util/common_util.py`
  - `cli/h_bonds/make_pp_matrix_ip6.py` (its binary parser)
- [x] B4. Update `make_pp_matrix_ip6.py` to consume canonical parser only.
- [x] B5. Add oxygen mode option to `make_pp_matrix_ip6.py` internals (`all`/`terminal`).
- [x] B6. Wire oxygen mode into CLI in `cli/app.py` for `pp-matrix`.
- [x] B7. Update `h_bonds_comparison_ip6_microstates.py` to use canonical parser path.
- [x] B8. Update `h_bonds_ip6.py` to use canonical parser path (remove local parser copies or wrap temporarily).

Acceptance criteria for Phase B:

- All H-bond code paths use one parser implementation.
- `IP_101101` recomputed matrix no longer shows implausible high `P4->P2`.
- Existing command behavior preserved except corrected values and new oxygen mode.

---

## Phase C - Controlled Cleanup (duplication and consistency)

- [ ] C1. Remove duplicated parser/analysis functions after migration:
  - `parse_xpm`, `analyze_hydrogen_bonds`, `parse_hbond_log_to_dataframe`, `load_hb_num_xvg`, `refine_atom_name` duplicates.
- [ ] C2. Consolidate phosphate mapping constants into one shared place for H-bond workflows.
- [ ] C3. Reconcile/centralize `draw_phosphorus_diagram` duplication (currently two implementations).
- [ ] C4. Ensure no remaining callers import removed duplicates.

Acceptance criteria for Phase C:

- `rg` confirms single authoritative implementation for each shared parsing function.
- Tests and CLI smoke runs pass.

---

## Phase D - Package-wide hygiene (non-blocking for bug fix)

- [ ] D1. Remove/archive old coordination summarizer module (`summerize_nap_dist_count_ip6_old.py`) from active paths.
- [ ] D2. Reduce `os.chdir` usage in CLI modules (prefer absolute paths).
- [ ] D3. Replace ad-hoc `print` statements with logging where practical.
- [ ] D4. Clarify CLI entrypoint docs (`cli.py` Typer stub vs `cli/app.py`).

Acceptance criteria for Phase D:

- No behavior regressions in command outputs.
- Cleaner dependency and logging behavior.

---

## 5) Execution Order (Recommended)

1. Phase A (tests first)
2. Phase B (correctness changes)
3. Recompute affected outputs for target project(s)
4. Phase C (remove duplicates)
5. Phase D (hygiene)

---

## 6) Verification Matrix

For each modified command:

1. `ip6-analysis hbonds ...`  
   - Verify generated `*_hb.log`, `*_hb_matrix.xpm`, downstream row-index consumers still align.

2. `ip6-analysis pp-matrix ...`  
   - Verify `*_PtoP_matrix.csv/.npy` values for known case(s).
   - Verify `species_*_union.pdf` and `species_*_sum.pdf` still render.
   - Verify `--oxygen-mode terminal` works and differs from `all`.

3. `ip6-analysis hbonds-compare ...`  
   - Verify summary values remain coherent after parser unification.

---

## 7) Data Rebuild Impact

Expected affected artifacts (must be regenerated where needed):

- `process/IP_*/analyze_*/h_bonds/*_PtoP_matrix.csv`
- `process/IP_*/analyze_*/h_bonds/*_PtoP_matrix.npy`
- `process/IP_*/analyze_*/h_bonds/species_*_union.pdf`
- `process/IP_*/analyze_*/h_bonds/species_*_sum.pdf` (if mode generated)
- Any tables/plots that use row-indexed `hbonds_per_index` via old parser paths

---

## 8) Rollback / Safety

1. Keep commits small and scoped by phase.
2. After each phase, run targeted checks before continuing.
3. If a phase fails verification, revert only that phase commit(s), not whole branch.

---

## 9) Work Log (update during execution)

Use this section as a live changelog.

### 2026-03-16

- Status: In progress (Phase A completed except one optional coverage item; Phase B completed through parser unification + oxygen mode).
- Files changed:
  - `src/electrofit_analysis/structure/util/hbond_io.py`
  - `src/electrofit_analysis/structure/util/common_util.py`
  - `src/electrofit_analysis/cli/h_bonds/h_bonds_ip6.py`
  - `src/electrofit_analysis/cli/h_bonds/make_pp_matrix_ip6.py`
  - `src/electrofit_analysis/cli/app.py`
  - `tests/conftest.py`
  - `tests/test_hbond_io.py`
  - `tests/test_pp_matrix_oxygen_mode.py`
- Commands run:
  - `python -m compileall src/electrofit_analysis`
  - `python -m pytest -q tests`
  - direct numerical validation scripts for `IP_101101` using patched parser/builder
- Validation results:
  - Test suite: `3 passed`
  - `IP_101101` corrected with aligned rows:
    - `P4->P2` (union, all): `0.00029997` (was `0.319068` in old artifact)
    - `P4->P2` (union, terminal): `0.00029997`
- Open issues:
  - Phase C cleanup still pending (full removal of duplicated wrappers and mapping centralization).

---

## 10) Current Status

- Plan created: 2026-03-16
- Overall state: `execution_in_progress`
- Next action: Phase C controlled cleanup (remove remaining wrappers/duplication) and regenerate affected project artifacts.
