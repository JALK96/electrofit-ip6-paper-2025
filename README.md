# electrofit-ip6-paper-2025

This repository contains the code to reproduce the
IP6 parameterization and subsequent analyses presented in the 2025 paper ....

- Primary code lives in: `electrofit-workspace/`
- See `envs/` for conda environments
- Usage instructions: see below

## Quickstart
Clone the repository,

```bash
git clone https://github.com/JALK96/electrofit-ip6-paper-2025.git
```
and install via mamba:
```bash
cd electrofit-workspace
mamba env create -f envs/electrofit-analysis.yml
conda activate electrofit-analysis
```
or using conda only:
```bash
cd electrofit-workspace
conda env create -f envs/electrofit-analysis.yml
conda activate electrofit-analysis
```

If you only plan on using the package to parameterize your own set of partial charges, you can install:
```bash
cd electrofit-workspace
conda env create -f envs/electrofit.yml
conda activate electrofit
```
If you prefer to install electrofit in an existing environment, you can do so by running:
```bash
cd electrofit-workspace
pip install -e ./packages/electrofit
```
and optionally
```bash 
pip install -e ./packages/ip6-analysis
```
to install the IP6 analysis package, to reproduce the analysis presented in the paper...
