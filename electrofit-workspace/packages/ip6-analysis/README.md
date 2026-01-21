## ip6-analysis (package notes)

This package provides the `ip6-analysis` CLI (`electrofit_analysis.cli.app`).

### `ip6-analysis coordination`

Runs an explicit coordination/count analysis (boolean tensor + per-phosphate counts)
and produces RDF plots for ion–peripheral oxygen distances.

#### Ion selection

Use `--ion-name` to select the ion by **atom name** (e.g. `NA`, `K`, `CA`, `MG`).

#### Project vs collection mode

Some datasets are organized as **multiple projects** (e.g. concentrations) inside one
ion-type folder:

```
K/
  50mM/   (project root; contains process/)
  100mM/
  150mM/
```

Use `--project-mode` to tell the CLI how to interpret `-p/--project`:

- `--project-mode project`: `-p` points to a single project root (expects `<project>/process/...`).
- `--project-mode collection`: `-p` points to a directory containing multiple projects
  (each child has its own `process/`).
- `--project-mode auto` (default): detect based on whether `process/` exists under `-p`.

#### RDF-derived cutoffs (default) and overrides

The explicit coordination/count algorithm needs a cutoff radius (nm).

- Default (`--coord-cutoff-nm` **not** provided): the cutoff is estimated from the RDF
  first-shell end (first minimum after the first peak), averaged across the six
  phosphates **within each microstate**.
- Fixed override: `--coord-cutoff-nm 0.32` forces a fixed cutoff (useful for reproducibility).
- Global cutoff across the scan scope: `--coord-cutoff-nm-global` computes one cutoff
  from all detected first-shell ends across the scan scope (project or collection)
  and applies it to all explicit counts in that run.

When global modes are enabled, it is recommended to also provide `-m/--molecule` to avoid
mixing unrelated microstates within the same `process/` tree.

#### RDF plots

RDF plotting options:

- `--rdf-show-first-shell`: draw a dotted guide at the detected first-shell end.
- `--rdf-cutoff-mode fixed|first-shell`: controls the RDF integration/guide used for the
  coordination number displayed in the RDF panels.

