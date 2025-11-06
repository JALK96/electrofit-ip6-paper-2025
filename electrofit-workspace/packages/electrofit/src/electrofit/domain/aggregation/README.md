# Step 6 — Ensemble Charge Aggregation and Symmetry Averaging

This module (`average_charges.py`) implements the “Step 6” aggregation phase of the
ElectroFit pipeline. It collects per‑conformer partial charges, computes ensemble
averages, optionally applies group (symmetry) averaging, updates a MOL2 with the
final charges, and runs ACPYPE to generate GROMACS topology.

The code is written to be deterministic and to keep orchestration logic thin;
plotting side effects are isolated into a dedicated `results/plots/` subfolder so
root‑level charge artifacts remain stable and unambiguous.

---

## Inputs and Configuration

The step operates inside a molecule process directory (`process/<mol>/`) and
expects the following structure and files:

- `extracted_conforms/` — per‑conformer outputs from Step 4/5 containing RESP‑fitted MOL2s.
- `run_gau_create_gmx_in/<mol>.acpype/` — folder with the input MOL2 (for atom names/type) and ACPYPE reference.
- `run_gau_create_gmx_in/*.json` or `extracted_conforms/equiv_groups.json` — optional symmetry groups mapping used for group averaging.

Effective configuration is resolved via the standard snapshot layering and read
with `load_config(...)`. Relevant fields are:

- `project.protocol` — `bcc` or `opt` (logged for context only).
- `project.charge` — target net charge (used for re‑normalisation after outlier removal).
- `project.atom_type` — e.g., `gaff2` (selects the source MOL2 naming).
- `symmetry.ensemble` — per‑stage mode, resolved to `(adjust_symmetry, ignore_symmetry)`
  for logging and plotting context; does not block group averaging.
- `project.calculate_group_average` — if `true`, apply symmetry group averaging when a symmetry JSON is present.
- CLI/runtime flags to this function: `remove_outlier`, `plot_histograms`, `hist_combine_groups`, `hist_bins`.

---

## Processing Stages (high‑level)

1. Compose a `results/` snapshot (`electrofit.toml`) under `process/<mol>/results/`.
2. Resolve configuration and symmetry flags for the ensemble stage (used for logging/plots).
3. Collect per‑atom charge series across conformers → `charges_dict.json` and atom‑name adjusted dictionary.
4. Compute the raw ensemble mean for each atom and persist:
   - `average_charges.chg` — one value per atom in MOL2 atom order.
   - `average_charges.txt` — human‑readable `#Atom_Name\tAverage_Charge` table.
5. Optional outlier removal (IQR‑based): recompute means on the cleaned set,
   re‑normalise to the target net charge, then persist:
   - `cleaned_adjusted_charges.json`
   - `cleaned_average_charges.chg`
   - `hist_no_outlier.pdf`
6. Locate a symmetry JSON (results → extracted_conforms → run directory).
7. If `calculate_group_average = true` and a mapping is available, apply group
   averaging to the base dataset (cleaned if outliers removed, else raw):
   - No outliers: `group_average_charges_dict.json`, `group_average_charges.txt`, `group_average_charges.chg`
   - With outliers removed: `cleaned_adjusted_group_average_charges_dict.json`,
     `cleaned_adjusted_group_average_charges.chg`
8. Update the input MOL2 using the best available `.chg` and run ACPYPE (user charge mode). See selection order below.
9. Emit plots under `results/plots/` (never clobber root charge files).
10. Log a decision summary (`aggregation` stage) including whether group averaging was applied.

---

## Final Charge Selection Order

Only one set of charges is embedded into the MOL2 and passed to ACPYPE. The
module picks the best candidate that exists in this priority order:

1. `cleaned_adjusted_group_average_charges.chg`
2. `cleaned_average_charges.chg`
3. `group_average_charges.chg`
4. `average_charges.chg`

This provides intuitive behaviour across combinations of outlier removal and
symmetry group averaging.

---

## Files Written (root of `process/<mol>/results/`)

Always present:
- `charges_dict.json` — per‑atom arrays of conformer charges with `average_charge` fields.
- `initial_charges_dict.json` — initial charges parsed from the source MOL2 for reference.
- `average_charges.chg` — raw ensemble mean (never overwritten by plots).
- `average_charges.txt` — the same raw means as a labelled table.

Outlier removal (optional):
- `cleaned_adjusted_charges.json`, `cleaned_average_charges.chg`, `hist_no_outlier.pdf`.

Group averaging (optional; requires symmetry JSON):
- `group_average_charges_dict.json`, `group_average_charges.txt`, `group_average_charges.chg`
  (or cleaned variants when outliers are removed).

MOL2/ACPYPE:
- `averaged_<mol>.mol2` (or `_cleaned` variant) — MOL2 updated with the selected `.chg` file.
- ACPYPE outputs are written next to this MOL2 (via `run_acpype(..., charges="user")`).

Manifest / Logs:
- `hist_manifest.json` (when histogram generation is enabled).
- `process.log` contains a decision summary and per‑file status.

---

## Plots (`process/<mol>/results/plots/`)

Plots are generated using the dictionaries described above; they never modify
root‑level `.chg` files:

- `charges.pdf` — violin plots of per‑atom distributions (raw vs. optional comparison).
- `charges_by_symmetry.pdf` — distributions coloured by symmetry groups.
- `charges_ip6.pdf` / `charges_comparison_ip6.pdf` — symmetry‑aware, IP6‑style atom plots.
- `hist.pdf` — per‑atom histograms before filtering.
- `hist_adjusted.pdf` — histograms with mean lines for the currently active charge set
  (ensemble averages or group averages).
- `hist_groups.pdf` — group‑combined histograms when `hist_combine_groups=true`.
- Plot helpers may also write `average_charges.chg` (and comparison variants) into this folder
  to record the mean values used in the figure; these are informational and do not affect the
  root‑level final charge selection.

---

## Decision Semantics

- `symmetry.ensemble` → `(adjust_symmetry, ignore_symmetry)` is logged for context
  (RESP stage behaviour), and controls which plots are produced; it does not block
  group averaging in Step 6.
- `calculate_group_average` → if `true` and a symmetry JSON is available, the module
  computes and uses group‑averaged charges (or cleaned group averages when outliers were removed).
- `remove_outlier` → applies IQR filtering before averaging and re‑normalises the net charge.

Corner cases:
- If raw ensemble means already satisfy symmetry relations, `group_average_charges.chg`
  will numerically match `average_charges.chg`.
- If `calculate_group_average=true` but no symmetry JSON is found, group averaging is skipped
  with a warning; ensemble averages remain the final charge set.

---

## Developer Notes

Key helpers:
- `_find_symmetry_json(...)` — best‑effort discovery of `equiv_groups.json`.
- `_apply_group_average_in_memory(...)` — in‑memory group averaging over any charge dict.
- `_write_charge_bundle(...)` — persists a complete `.json/.txt/.chg` triple for a charge dict.

This module aims to make the dataflow explicit and easy to reason about:
- Root‑level `.chg` files are authoritative candidates for the final MOL2 update.
- Plotting outputs are sandboxed to `results/plots/` and do not influence the final choice.
- One ACPYPE run is executed per aggregation call, using the best available `.chg` per the selection order.

If you need to add another averaging or filtering mode, prefer building a new
charge dict and persisting it via `_write_charge_bundle(...)`, then plug it into the
final selection list with a clear priority.

