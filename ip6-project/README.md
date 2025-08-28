# IP6 project

This directory is the project root for the IP6 microstate parameterization
and analysis. It contains the inputs used for parameterization, the per‑
microstate pipeline outputs and the post-processed/aggregated results.

## Top-level layout

- `data/` — input materials and auxiliary files used by the pipeline:
  - `data/input/<MICROSTATE>/` — per-microstate input files (e.g. `<name>.mol2`,
    `<name>.pdb`, `electrofit.toml`, `equiv_groups.json`, `input.ef`).
  - `data/MDP/` — model MDP templates used for GROMACS runs.
- `process/` — per-microstate runtime outputs and results. Each microstate has
  its own directory (e.g. `process/IP_111101/`).

There may also be auxiliary files such as `README.txt` in `data/` containing
project-specific notes.

## Per-microstate `process/<MICROSTATE>/` layout

Each microstate directory contains subfolders that reflect the
pipeline steps and post-processing. Example structure for `process/IP_111101/`:

- `run_gau_create_gmx_in/` — initial Gaussian / antechamber / acpype inputs and
  generated artefacts; contains `*.mol2`, `*.acpype/` directories and
  a `process.log` for this stage.
- `run_gmx_simulation/` — GROMACS setup and simulation outputs (topologies,
  trajectories, logs).
- `extracted_conforms/` — extracted conformers; contains many conformer
  subdirectories `IP_111101c0/`, `IP_111101c1/`, ... Each conformer folder
  usually contains:
  - `IP_...cN.pdb` — conformer coordinates
  - `IP_...cN.gcrt`, `IP_...cN.gesp`, `IP_...cN.gcrt.log` — Gaussian inputs/outputs
  - `IP_...cN.mol2`, `IP_...cN_resp.mol2` — mol2 before/after RESP
  - `IP_...cN-resp1.*`, `IP_...cN-resp2.*` — RESP stage outputs
  - `process.log` — stage log for that conformer
  - `symmetry_resp_MOD.txt` / `equiv_groups.json` when symmetry was applied
- `results/` — aggregated and final per-microstate results of the parameterization. Typical files:
  - `averaged_<MICROSTATE>.mol2`, `average_charges.chg`, `group_average_charges.chg`
  - `charges.pdf`, `charges_by_symmetry.pdf`, `plots/` — visual summaries
  - `charges_dict.json`, `group_average_charges_dict.json` — machine-readable summaries
  - `results.log` — log of the aggregation/update step
- `run_final_gmx_simulation/` and `analyze_final_sim/` — final simulations
  (with parameterized charges) and their analysis outputs.

## How to rerun the parameterization / analysis

To rerun the parameterization or analysis for a single microstate, make a copy
of the input and run the pipeline steps.

1. Prepare a fresh project directory

    ```bash
    mkdir -p /path/to/ip6-rerun/data/input/
    ```

2. Copy the microstate inputs and MDP templates

    ```bash
    # copy one microstate (replace <MICROSTATE> and adjust paths as needed)
    cp -a ip6-project/data/input/<MICROSTATE> /path/to/ip6-rerun/data/input/

    # copy MDP templates used for simulations
    cp -a ip6-project/data/MDP /path/to/ip6-rerun/data/
    ```

3. Run the electrofit pipeline steps
  (recommended: run each step and inspect logs before continuing)

    ```bash
    electrofit step0 --project /path/to/ip6-rerun
    electrofit step1 --project /path/to/ip6-rerun
    electrofit step2 --project /path/to/ip6-rerun
    electrofit step3 --project /path/to/ip6-rerun
    electrofit step4 --project /path/to/ip6-rerun --sample 100
    electrofit step5 --project /path/to/ip6-rerun
    electrofit step6 --project /path/to/ip6-rerun
    ```

    Tip:
    - Run with verbose logging to capture more output (helpful for debugging), e.g.:
  
    ```bash
    ELECTROFIT_LOG_LEVEL=DEBUG \
    electrofit step5 --project /path/to/ip6-rerun
    ```

4. Locate the generated partial charges and aggregated results

    After the pipeline completes, per-microstate results live under
    `/path/to/ip6-rerun/process/<MICROSTATE>/results/` (for example
    `averaged_<MICROSTATE>.mol2`, `average_charges.chg`, `results.log`).

5. (Optional) Final simulations

    If you want to run the final GROMACS simulation with the new charges,
    execute:

    ```bash
    electrofit step7 --project /path/to/ip6-rerun
    electrofit step8 --project /path/to/ip6-rerun
    ```

6. (Optional) Post-processing with `ip6-analysis`

    If `ip6-analysis` is installed in the same environment as `electrofit`, you
    can run post-processing commands such as:

    ```bash
    ip6-analysis hbonds --project /path/to/ip6-rerun
    ip6-analysis pp-matrix --project /path/to/ip6-rerun
    ip6-analysis coordination --project /path/to/ip6-rerun
    ```

7. Notes:

    - If you want to rerun the entire set of microstates, you can copy the entire input directory from `ip6-project/data` to `/path/to/ip6-rerun/data`

    - After copying, you can run the pipeline steps as described above. The steps will process all microstates in the new input directory, without specifying individual microstates.
  
    - You would still be able to run steps for individual microstates by using the `--molecule` flag, e.g.

    ```bash
    electrofit stepX --project /path/to/ip6-rerun --molecule <MICROSTATE>
    ```

    - the `<MICROSTATE>` placeholder should be replaced with the actual microstate name (e.g. `IP_010101`).
