# electrofit (core)

Minimal core tooling to process small molecules / residues, run short MD (GROMACS) and extract conformers for charge workflows.

## Quick Start

Create / activate an environment, install the package (editable example):

```bash
pip install -e ./packages/electrofit
```

Prepare a project folder (contains `electrofit.toml`, `data/input/<MOLNAME>/<files>`):

```bash
electrofit step0 --project /path/to/project
```

Run successive workflow steps (see below). Each step reads the project root you pass with `--project` and optional config overrides via `--config`.

## Configuration (`electrofit.toml`)

Key sections (example):

```toml
[project]
molecule_name = "IP_011101"
residue_name  = "IP6"
charge        = -8
protocol      = "bcc"     # or "opt"

[paths]
mdp_dir = "data/MDP"
base_scratch_dir = "/scratch/${USER}/electrofit"

[gromacs.runtime]
threads = 8
pin = true

[sampling]
method = "maxmin"   # linear | random | maxmin
count = 25
seed = 123
```

CLI arguments override config values. (Früherer ENV-Fallback `ELECTROFIT_CONFIG_PATH` wurde entfernt – nur noch explizit via `--config`.)

## Workflow Steps (0–4)

### step0 – Initial project scaffold
Creates baseline directory structure (e.g. `process/`), validates inputs, copies / normalises initial structural files if required.

Command:
```
electrofit step0 --project <project_root>
```

### step1 – Initial processing
Legacy pre-processing (Gaussian / RESP input staging, symmetry handling). Produces intermediate files placed under `process/<mol>/run_gau_create_gmx_in`.

### step2 – Simulation directory setup
Builds per-molecule GROMACS run directories (`run_gmx_simulation`), inserts MDP templates, generates topology / coordinate starting point.

### step3 – Start short MD production
Runs energy minimisation + equilibration + (short) production with parameterised runtime flags (threads, pinning). Produces `md_center.xtc` and `md.gro` used later.

```bash
electrofit step3 --project <project_root>
```

### step4 – Extract representative conformers
Reads each `process/<mol>/run_gmx_simulation` trajectory and writes sampled conformers into `process/<mol>/extracted_conforms/` (PDB + copied ancillary input files).

```bash
electrofit step4 --project <project_root> \
	--sample 20 --sampling-method maxmin --seed 123 --workers 2
```

## Conformer Sampling (Step 4)

Supported strategies (`--sampling-method`):

* `linear` (default): evenly spaced frames (deterministic)
* `random`: uniform random without replacement (seeded)
* `maxmin`: farthest-point (max–min) diversity in RMSD space

Important options:

* `--sample N` desired number of conformers (fallback 20 or `[sampling].count`)
* `--sampling-method METHOD` selection strategy
* `--seed S` seed for random / starting frame in maxmin
* `--workers K` parallelise across molecule directories (0 = auto heuristic)
* `--no-progress` disable the progress bar

Config section (optional):

```toml
[sampling]
method = "maxmin"
count = 25
seed = 123
```

CLI overrides config. Output line includes method, e.g. `Extracted 25 conformers (method=maxmin) ...`.

Performance notes:

* `maxmin` costs roughly O(k * N) RMSD evaluations (k = selected frames). For very large trajectories try a smaller `--sample` first.
* Parallelisation is at the granularity of molecules, not frames inside one trajectory.

Planned extensions: RMSD cutoff incremental selection, k-medoids clustering, PCA/grid stratified sampling.

## Parallel Execution

`--workers` uses a process pool. If unset / 0 an automatic value (CPU count minus one) capped by number of molecules is used. Sequential mode (workers=1) retains a live progress bar.

## Logging & Scratch

Scratch directories are created under `base_scratch_dir` (environment variables expanded). Finalisation is idempotent. GROMACS step automates interactive selections (e.g. `trjconv`).

## Architecture Overview & Ongoing Refactor

Layered Struktur (in Arbeit):

```
electrofit/
	domain/              # Fachlogik (prep, charges, symmetry, sampling)
	adapters/            # Technische Tool-Anbindungen (gromacs, gaussian, resp)
	infra/               # Infrastruktur (scratch, config snapshot, logging)
	workflows/           # Thin CLI-Orchestrierung (step*-Skripte)
	io/                  # Dateiformat Operationen (mol2, resp, files)
	viz/                 # Plotting / Visualisierung (optional)
```

Zentrale Prinzipien:
* Domain-Funktionen kapseln Logik ohne direkte Shell-Kommandos.
* Adapters kapseln externe Tools (GROMACS bereits migriert; Gaussian/RESP folgen).
* Einheitliche Deprecation-Shims für alte `core.*` Entry Points (einmalige Warnung).

### Gemeinsame RESP Symmetry Utility

Vorher doppelte Logik in Prep & Charges: Jetzt zentralisiert in `domain.symmetry.resp_constraints.apply_and_optionally_modify`.
Vorteile: Single Source of Truth für JSON Equivalence Groups Anwendung; konsistente Logging & Fehlerbehandlung.

### Geplanter Adapter-Ausbau

Neu (Skeletons):
* `adapters/gaussian.py` – Bau & Ausführung von Gaussian Inputs + Cache-Hydration API.
* `adapters/resp.py` – Zwei-Phasen RESP Pipeline (ac/respgen + symmetry + charge update).

Aktuelle Domain-Pfade rufen noch `cli.run_commands` Funktionen direkt; schrittweise Migration wird diese durch Adapter-Funktionen ersetzen, wodurch Testbarkeit (Mocking / Timing) verbessert und zukünftige Engines (z.B. Psi4) pluggable werden.

Migrationsstrategie (Kurz):
1. Skeletons hinzufügen (done).
2. Domain `process_conformer` & `process_initial` intern umstellen auf Adapter-Aufrufe (keine Public API Änderung).
3. Entfernen direkter `run_command`-Abhängigkeiten aus Domain (nur Adapters dürfen shellen).
4. Optional: strukturierte Rückgabeobjekte (Timing, Pfade, Checksums) für reproducibility / Caching.

Siehe `docs/REFACTOR_PLAN.md` für detaillierten Status & weitere Phasen.

## Diagnostic Histograms & Outlier Filtering (Step 6)

Optional diagnostic visuals and IQR-based outlier removal for charge ensembles.

Enable via flags (independent / combinable):

```
electrofit step6 --project <proj> \
	--plot-histograms \
	--remove-outlier \
	--hist-combine-groups \
	--hist-bins 30 \
	--outlier-iqr-factor 1.5
```

Artifacts (within `process/<mol>/results/`):

| File | Description |
|------|-------------|
| hist.pdf | Per-atom distributions before filtering |
| hist_no_outlier.pdf | Overlay before vs after outlier removal |
| hist_adjusted.pdf | After applying group/symmetry average (or cleaned reweighted) |
| hist_groups.pdf | Group-combined distributions (symmetry groups merged) |
| hist_summary.json | Structured report (removed counts, per_atom_outliers, indices, IQR factor) |
| hist_manifest.json | Machine-readable status for each optional histogram (expected/created/reason) |
| cleaned_adjusted_charges.json | Filtered charges per atom |
| cleaned_adjusted_group_average_charges* | Group-averaged filtered charges |

Outlier criterion (per atom column): value < Q1 - f*IQR or > Q3 + f*IQR (factor f = `--outlier-iqr-factor`, default 1.5). A conformer is discarded if any atom is an outlier (union mask) – conservative pruning to remove multi-atom anomalies.

Sequence (if all enabled):
1. hist.pdf (raw)
2. hist_no_outlier.pdf (after removal)
3. hist_adjusted.pdf (after group averaging / reweight)
4. hist_groups.pdf (aggregated symmetry groups)

Design rationale:
* Zero-impact unless flags are passed.
* JSON summary supports automated regression thresholds.
* Separation of diagnostics from core pipeline ensures reproducibility.

Planned extensions: z-score alternative, adaptive IQR, persistent whitelist of conformers.

### Histogram Manifest (`hist_manifest.json`)

Because some histogram artefacts are *conditional* (e.g. `hist_adjusted.pdf` only makes sense once group averaging or adjusted means are applied), the pipeline writes a manifest capturing — per logical histogram stage — whether it was expected and whether it was actually created.

Schema (per key):

```jsonc
{
	"initial": { "expected": true,  "created": true,  "path": "hist.pdf",              "reason": null },
	"after_outlier": { "expected": true,  "created": true,  "path": "hist_no_outlier.pdf", "reason": null },
	"adjusted": { "expected": true,  "created": false, "path": null, "reason": "no groups" },
	"group_combined": { "expected": false, "created": false, "path": null, "reason": "hist_combine_groups flag not set" }
}
```

Interpretation rules:
* `expected=false` means the user did not request (via flags) or preconditions were deliberately not met (e.g. symmetry disabled) — missing artefact is **not** an error.
* `expected=true` & `created=true` implies a file at `path` exists.
* `expected=true` & `created=false` signals a soft failure or structural reason; `reason` provides diagnostics (e.g. `no non-empty series`, `error: <msg>`, `no groups`). Tests may choose to fail in CI if such cases become frequent.
* Stable contract: new histogram types may be appended (additive change) — consumers should ignore unknown keys for forward compatibility.

Rationale: Avoid brittle test expectations for artefacts that depend on small synthetic datasets or rare branches while still surfacing why something was skipped. This pattern is preferable to writing placeholder PDFs that convey no information.

Testing Guidance:
* Unit tests assert required base artefacts (`hist.pdf`, `hist_no_outlier.pdf`) exist when requested.
* Conditional artefacts (`hist_adjusted.pdf`, `hist_groups.pdf`) are validated through the manifest (`expected -> created`).
* When adding new histogram stages, update the manifest keys and extend tests only if the new stage should always be produced under current flags.

Backward Compatibility: Legacy scripts looking only for `hist.pdf` keep working; the manifest is an additive enhancement.

## Testing

Unit tests include sampling selection determinism and basic functional checks. Shortened MDP templates keep test runtime low.

## License

## Konfig-Präzedenz & Fill-In

Die endgültige `electrofit.toml` Snapshot-Datei in jedem Arbeitsverzeichnis entsteht über einen klar definierten Layering-Prozess:

Starke Overrides (überschreiben vorhandene Werte in Reihenfolge – spätere gewinnt):
1. Molekül-spezifische Eingabe: `data/input/<MOL>/electrofit.toml`
2. Prozess-spezifische Datei: `process/<MOL>/electrofit.toml` bzw. vorherige Schritt-Ausgabe (z.B. `run_gau_create_gmx_in/electrofit.toml` oder `results/electrofit.toml` je nach Schritt)
3. CLI `--config` (falls angegeben)

Fill-In Ebene (füllt nur fehlende Schlüssel, überschreibt niemals):
4. Projektweite Defaults: `<project_root>/electrofit.toml`

Semantik:
* „Override“ ersetzt bestehende Werte (deep merge, scalars & Subtrees vollständig überschrieben).
* „Fill-In“ ergänzt nur Keys, die noch nicht existieren (rekursiv), lässt vorhandene Werte unverändert.
* Dieses Muster verhindert, dass projektweite Defaults unabsichtlich molekül-spezifische Einstellungen verdrängen, reduziert aber Duplikation bei globalen Parametern.

Beispiel – gewünschte Kraftfeld-Priorität:
* Molekül: `simulation.forcefield = "amber14sb.ff"`
* Projekt setzt keinen Forcefield-Schlüssel → Snapshot übernimmt Molekülwert.
* Falls Projekt später `simulation.forcefield` hätte, würde er NICHT das molekül-spezifische überschreiben (weil Projekt nur Fill-In ist) – gewünschtes Verhalten.

Logging-Markierungen:
* `[config][override] key: old -> new` für Überschreibungen durch starke Layer.
* `[config][fill] key: value` wenn ein fehlender Schlüssel durch Projekt-Defaults ergänzt wurde.
* `[config] no overrides applied ...` / `no fills` wenn keine Änderung.

CLI `--config`:
* Wird als stärkste Override-Ebene eingehängt (nach Molekül & Prozess), um gezielt Werte temporär zu ersetzen.
* Überschreibt vorhandene Werte, führt keine Fill-Ins durch.

Implementierung:
* Zentral in `compose_snapshot` (`infra/config_snapshot.py`, legacy Alias `workflows/snapshot.py`).
* Alle relevanten Schritte (1,2,3,4,6,7) benutzen diese Funktion, um doppelte Logik zu vermeiden.

Vorteile:
* Vorhersehbare Priorität ohne nachträgliche „Hack“-Korrekturen.
* Minimierte Redundanz in per-molekül TOMLs (nur spezifisches definieren, globales bleibt im Projektfile).
* Klare Trennung zwischen „Ändern“ (Override) und „Auffüllen“ (Fill-In).

Best Practices:
* Molekül-spezifische Parameter (Ladung, Forcefield, besondere Sampling-Optionen) in der Molekül-Datei definieren.
* Einheitliche globale Pfade / Laufzeit-Defaults nur im Projektroot pflegen.
* `--config` für experimentelle Runs oder CI-spezifische Tuning verwenden (nicht dauerhaft einchecken).

Edge Cases:
* Fehlt ein Molekül-TOML komplett, greifen Prozess/CLI/Projekt in genannter Reihenfolge.
* Mehrere Moleküle: bestimmte restriktive Keys werden geloggt, wenn sie in Multi-Mol Kontext überschrieben würden (siehe `RESTRICTED_DEFAULT`).

TBD (add your license information here).
