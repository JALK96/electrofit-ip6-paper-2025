# ElectroFit

ElectroFit derives partial charges for flexible molecules from a conformational
ensemble. It samples conformers by MD, fits RESP charges to each conformer with
Gaussian, averages them (optionally over symmetry-equivalent atoms) and sets up
GROMACS simulations with the averaged charges. `ip6-analysis` analyses
simulations of myo-inositol hexakisphosphate (IP6): ion coordination and
distances, hydrogen bonds, charges and dihedrals.

## Contents

- `electrofit-workspace/packages/electrofit/`: the ElectroFit package and its
  `electrofit` command ([README](electrofit-workspace/packages/electrofit/README.md)).
- `electrofit-workspace/packages/ip6-analysis/`: the IP6 analysis package and its
  `ip6-analysis` command ([README](electrofit-workspace/packages/ip6-analysis/README.md)).
- `electrofit-workspace/data/`: force fields (`amber14sb.ff`,
  `amber14sb_corrected-Na.ff`), MDP templates and example inputs.
- `electrofit-workspace/examples/`: example projects.
- `electrofit-workspace/tests/`: tests.
- `electrofit-workspace/envs/`: conda environment definitions.
- `ip6-project/data/`: inputs and MDP templates for the 11 IP6 protonation
  microstates at pH 7.4 ([README](ip6-project/README.md)).
- `notebooks/traj.ipynb`: trajectory visualisation; set the data paths in its
  first cell. Its environment is defined in `traj_vis_env.yml`.

## Installation

From `electrofit-workspace/`:

```bash
conda env create --prefix /path/to/env --file envs/electrofit-dev.lock.yml
conda activate /path/to/env
python -m pip install --no-deps -e ./packages/electrofit -e ./packages/ip6-analysis
electrofit --help
ip6-analysis --help
```

Gaussian, AmberTools (antechamber, ACPYPE) and GROMACS are external programs.
Check that the versions available on the host work with the workflow before a
full parameterization or MD run.

## Running

Each `electrofit` step takes a project directory (`--project`); see the
[ElectroFit README](electrofit-workspace/packages/electrofit/README.md) for the
steps and configuration. `electrofit-workspace/scripts/run/run_all_*.sh` run the
steps for one molecule of the integration-test project; `PROJECT` and `CACHE`
select another project and cache, and `DRY_RUN=1` prints the commands without
running them.
