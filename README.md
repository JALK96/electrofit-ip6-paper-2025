# ElectroFit software and starting inputs

Reusable parameterization and analysis source, with IP6 starting inputs.
The repository name is retained for now; it does not identify a published paper.

- `electrofit-workspace/packages/electrofit/`: parameterization implementation.
- `electrofit-workspace/packages/ip6-analysis/`: analysis implementation.
- `electrofit-workspace/{scripts,envs,docs,tests,data,examples}/`: supporting code,
  environment definitions, documentation, tests and reusable input/example files.
- `ip6-project/data/`: microstate inputs and MDP templates to copy into a fresh
  study run. See [parameterization steps](ip6-project/README.md).
- `notebooks/`: existing trajectory exploration notebook; its historical data
  paths need selecting or repairing for each study. `traj_vis_env.yml` describes
  its environment.

Keep simulation outputs in the consuming study, outside this source checkout.
Historical outputs and manuscripts remain in their original study/publication
locations and are not part of this working tree.

## Working environment

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate /home/johannal96/PhD/software/envs/electrofit-dev
electrofit --help
ip6-analysis --help
```

This development environment uses editable packages from this checkout. The
existing `electrofit-analysis` environment also uses these packages. The legacy
`electrofit` environment uses only ElectroFit from this checkout. Source edits
apply to subsequent commands in all environments using the edited package. Keep recorded production builds
separate. Gaussian, AmberTools/ACPYPE and GROMACS compatibility must be checked
for the intended host before running a full parameterization or MD workflow.

Dependency versions are recorded in
`electrofit-workspace/envs/electrofit-dev.lock.yml`. To create a separate working
environment, run from `electrofit-workspace/`:

```bash
conda env create --prefix /path/to/new-env --file envs/electrofit-dev.lock.yml
conda activate /path/to/new-env
python -m pip install --no-deps -e ./packages/electrofit -e ./packages/ip6-analysis
```

The sparse checkout includes software, notebooks and starting inputs. Inspect
its selection with `git sparse-checkout list`. Keep simulation outputs in the
consuming study; `scripts/run/run_all_*.sh` accept `PROJECT` and `CACHE` overrides.
Use `DRY_RUN=1` to inspect those command sequences.
