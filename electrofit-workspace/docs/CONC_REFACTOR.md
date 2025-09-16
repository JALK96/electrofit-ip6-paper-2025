# Ion Concentration Refactor Plan

## Background
- Current configuration exposes `[simulation.ions]` with keys `cation`, `anion`, and
  `concentration`. The value is handed directly to `gmx genion -conc`, after an
  `editconf -d` expansion. The workflow therefore assumes:
  - The box edge is controlled via `editconf -d`, i.e. padding in nm around the solute.
  - Ion counts are whatever `genion` computes for the instantaneous box volume plus
    the extra ions required by `-neutral`.
- Users now need two modes:
  1. **Salt concentration mode** (legacy) – specify the molarity of salt pairs, keep
     the existing pipeline intact.
  2. **Fixed ion count + concentration mode** – request a specific number of ions of
     given species (e.g. 12 Na⁺) *and* a target concentration (e.g. 150 mM). To
     satisfy both, the simulation cell has to be sized so that the requested count is
     consistent with the requested molarity. This must work for any supported box
     type (`cubic`, `triclinic`, `dodecahedron`, …).

## Desired UX
```toml
[simulation.ions]
cation              = "NA"
anion               = "CL"
salt_concentration  = 0.15   # mol/L (optional; legacy behaviour)
enforce_neutrality  = true

[simulation.ions.targets]
# optional advanced mode
[simulation.ions.targets.positive_ion]
desired_count  = 12
concentration  = 0.150  # mol/L for that species

[simulation.ions.targets.negative_ion]
desired_count  = 12
# concentration must be omitted here; asymmetric salt support is not available yet
```
- If `targets` is absent → do exactly what we do today, but the key is now
  `salt_concentration` (the loader should still honour the historic `concentration`
  key with a warning/deprecation shim).
- If targets are provided → compute the box volume that reaches the given molarity
  with the desired count(s), then build a box of that volume *before* solvation.
  `salt_concentration` and the legacy padding value `box.edge_nm` should both be
  ignored (with warnings) whenever targets are present; users must pick one mode or
  the other.
- Neutrality rule (when `enforce_neutrality = true`):
  - If only one ion type is specified, infer the counter-ion count as
    `N_missing = N_known - q_solute` (sign convention chosen so the total charge is
    zero). Validate `N_missing >= 0`.
  - If both ion types are specified, check that `N_pos - N_neg + q_solute = 0`. If the
    identity fails, abort with an actionable error message explaining the mismatch and
    reporting the corrected counts to place in the TOML.
  - When `enforce_neutrality = false`, keep the user-specified counts even if the
    system charge is non-zero, but emit a high-visibility warning.

## Approach Overview
1. **Config schema & compatibility**
   - Extend `SimulationIonsSection` with the new keys (`salt_concentration`,
     `targets`, `enforce_neutrality`).
   - Update the loader/parser to convert old configurations (`concentration`) into
     the new structure and emit a warning. Enforce mutual exclusivity between
     `salt_concentration`/`targets` and between `box.edge_nm`/target mode to avoid
     ambiguous setups.
   - Update `RESTRICTED_DEFAULT` to treat the new keys as restricted in multi-molecule
     merges.
2. **Internal representation**
   - Introduce a helper (e.g. `electrofit.simulation.ions` module) that produces a
     canonical description from the config: either `SaltConfig` or
     `IonTargetConfig`. This should encapsulate the calculations (volume, counts,
     neutrality checks) so the adapters remain thin.
3. **Neutrality and ion bookkeeping**
   - The helper determines final ion counts by combining solute charge, target
     entries, and `enforce_neutrality`:
       * Compute species individually if both positive and negative targets are
         provided; otherwise infer the missing species.
       * Validate `N >= 0` for every species. When `enforce_neutrality` is true,
         refuse setups that cannot satisfy neutrality (e.g. `N_fixed < |q_solute|`).
       * Surface warnings when neutrality is disabled or when counts are inferred.
       * When both species specify concentrations, raise a configuration error during
         load with instructions to remove one of the concentration entries; asymmetric
         salt handling will be implemented later.
4. **Box sizing logic**
   - Write a utility to compute edge lengths from a target volume for each box type:
       * `cubic`: `L = V^(1/3)`.
       * `triclinic`/`rectangular`: simple decomposition.
       * `dodecahedron`: use GROMACS’ rhombic dodecahedron volume relation
         (`V = (16/9) * sqrt(3) * a³`, where `a` is the short box vector). Confirm
         the formula via documentation or experiments.
   - For target mode, base the volume calculation on whichever species provides a
     concentration entry (currently only the positive ion). If neither carries a
     concentration, raise an error explaining that volume cannot be determined.
   - Provide the inverse helper as well: given edge length(s), report volume for
     logging and validation.
5. **GROMACS adapter changes**
   - Refactor `set_up_production` so geometry/ion logic lives near the top:
       * Accept a new dataclass/struct describing the ions mode (includes counts for
         both species, computed/inferred flags, target concentrations, etc.).
       * If in salt mode → behave as today, using `editconf -d` and `genion -conc`.
       * If in target-count mode → compute the desired vectors, emit
         `editconf -f m.gro -o ... -box ax ay az -angles alpha beta gamma` (if
         necessary), or call `gmx editconf -bt ... -box ...` depending on type.
       * In target mode, ignore the legacy padding parameter (`-d`) and `-conc`, and
         instead call `genion` with explicit `-np`/`-nn` derived from the helper.
       * Keep `-neutral` optional: when `enforce_neutrality` is false, allow the user
         to request non-neutral setups by skipping the flag.
   - Add logging to record: requested vs inferred counts, computed volume, resulting
     box vectors, and ions inserted (parse stdout from `genion`). Surface warnings when
     neutrality constraints are disabled or violated, or when `salt_concentration` or
     `edge_nm` were ignored because targets were supplied.
6. **Pipeline touch points**
   - `pipeline/steps/step3.py` and `domain/final_sim.py` currently extract
     `conc`/`d` directly. Update them to pass the structured config (including the
     neutrality flag) instead of individual floats. Ensure they respect the
     one-mode-at-a-time rule and issue warnings when a deprecated knob is ignored.
   - Ensure `compose_snapshot` propagates the new key names and defaults are updated
     in project scaffolding (`data/MDP`, sample `electrofit.toml`).
7. **Documentation & migrations**
   - Update docs/samples to show both configuration modes and the neutrality flag.
   - Outline migration guidance (e.g. rename `concentration` → `salt_concentration`).
   - Provide a warning in logs when legacy keys are read, when neutrality is disabled,
     or when `salt_concentration`/`edge_nm` are ignored due to targets.
8. **Testing & validation**
   - Add unit tests for config parsing, neutrality checks, and geometry helpers.
   - Add regression tests covering:
       * Fixed-count scenarios with positive, negative, and neutral solute charges.
       * Both species specified vs single-species inference.
       * Neutrality disabled, ensuring warnings are emitted but processing continues.
       * Configurations that attempt to set concentration on both species, verifying
         that we raise the intended error message.
       * Configurations that set both `salt_concentration`/`edge_nm` and `targets`,
         verifying we warn and honour targets.
   - Add integration coverage (dry-run) to assert that generated `editconf` commands
     match expectations for cubic and dodeca boxes, and that ion counts sent to
     `genion` reflect the helper output.

## Implementation Steps (Detailed)
1. **Schema update**
   - Modify `SimulationIonsSection` in `config/loader.py` to include
     `salt_concentration`, `targets`, and `enforce_neutrality`. Define nested
     dataclasses for `IonTarget` (desired_count, optional concentration flag,
     maybe charge sign).
   - Adjust `_merge_into_dataclass`/`dump_config` to handle the nested structures.
   - Update `config/legacy.py` so legacy `.ef` imports map `IonConcentration` to the
     new field and set `enforce_neutrality = true` by default.
2. **Loader compatibility**
   - After reading TOML, detect a plain `concentration` key; move it to
     `salt_concentration` and log a deprecation warning.
   - If both positive and negative ion targets specify `concentration`, raise a
     descriptive error instructing the user to remove one of them (asymmetric support
     pending).
   - If both `salt_concentration`/`edge_nm` and `targets` are provided, issue a warning
     that the legacy scalar controls will be ignored and proceed with target mode.
   - Update `RESTRICTED_DEFAULT` in `config/merger.py` to reference
     `simulation.ions.salt_concentration`, `simulation.ions.targets.*`, and
     `simulation.ions.enforce_neutrality`.
3. **Ion config helper**
   - Create `electrofit/simulation/ions.py` with functions:
       * `build_ion_config(cfg.simulation.ions, solute_charge) -> SaltMode | TargetMode`.
       * `compute_needed_volume(target_count, molarity) -> float`.
       * `compute_box_vectors(volume, box_type) -> tuple[vectors]`.
       * `infer_counter_ions(targets, solute_charge, enforce_neutrality)`.
   - Include validations (non-negative counts, molarity > 0, no neutrality conflicts).
   - When neutrality enforcement fails, raise an exception whose message includes the
     observed counts, solute charge, and suggested corrected counts for the TOML.
   - Return a structure that clearly indicates which counts were user-provided vs
     inferred to assist logging.
4. **Adapter refactor**
   - Update `set_up_production` signature to accept the new ion spec; keep backward
     compatibility path for `conc`/`d` until all callers are migrated.
   - Insert geometry calculation before calling `editconf`; log computed box
     dimensions.
   - In target mode: call `genion` with `-np`/`-nn` derived from the helper and omit
     `-conc` entirely. When neutrality is disabled, skip `-neutral` to prevent
     `genion` from altering counts.
5. **Call site changes**
   - Update `pipeline/steps/step3.py` and `domain/final_sim.py` to gather the new
     fields and pass them to `set_up_production` (likely via the helper).
   - Ensure CLI defaults (`templates/electrofit.toml`) expose the neutrality flag and
     illustrate both the simple and target-count configurations.
6. **Docs & messaging**
   - Update `docs/README.md` and other tutorials to explain configuration options and
     neutrality logic (including error/warning scenarios).
   - Mention the change in release notes / changelog and highlight new log outputs.
7. **Testing strategy**
   - Unit tests for config parsing, neutrality checks, and geometry conversions.
   - Regression tests verifying computed box edge for a cubic case (≈5.08 nm for 12
     ions @ 150 mM) and that counter-ion counts neutralise systems with positive and
     negative solute charges.
   - Tests ensuring warnings/errors trigger correctly when neutrality requirements
     are unmet or disabled, when dual concentrations are specified, and when both
     salt and target modes are set simultaneously.

## Open Questions & Follow-ups
- Confirm the exact formula GROMACS uses for `-bt dodecahedron` volumes.
- Decide how to handle asymmetric salts once we introduce per-species concentration
  support (current plan forbids it and raises on load).
- Determine whether neutrality violations should always abort or whether we offer an
  additional override flag to degrade errors to warnings (current plan favours
  strict enforcement unless `enforce_neutrality = false`).
- Consider persisting the computed box vectors and final ion counts in the snapshot
  TOML for auditability (optional but could help debugging).
