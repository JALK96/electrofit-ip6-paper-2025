# packages/electrofit/src/electrofit/cli/app.py
import argparse
import os
import runpy
import sys
import logging
from pathlib import Path

STEP_MODULES = {
    "step0": "electrofit.pipeline.steps.step0",
    "step1": "electrofit.pipeline.steps.step1",
    "step2": "electrofit.pipeline.steps.step2",
    "step3": "electrofit.pipeline.steps.step3",
    "step4": "electrofit.pipeline.steps.step4",
    "step5": "electrofit.pipeline.steps.step5",
    "step6": "electrofit.pipeline.steps.step6",
    "step7": "electrofit.pipeline.steps.step7",
    "step8": "electrofit.pipeline.steps.step8",
}

def _run_module(module: str, project: Path, config: Path | None, rest: list[str]):
    # Export project path only (old ELECTROFIT_CONFIG_PATH removed)
    os.environ["ELECTROFIT_PROJECT_PATH"] = str(project)

    prev_argv = sys.argv
    # Only pass --project; keep additional arguments
    sys.argv = [module, "--project", str(project)] + rest
    try:
        runpy.run_module(module, run_name="__main__")
    finally:
        sys.argv = prev_argv

def _ensure_project_scaffold(project: str | os.PathLike) -> tuple[Path, Path]:
    root = Path(project).resolve()
    process = root / "process"
    if not process.exists():
        logging.info(f"'process' folder not found at {process}; creating it for you.")
        process.mkdir(parents=True, exist_ok=True)
    return root, process

def main():
    # Be talkative by default unless caller configured logging already
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(message)s")

    p = argparse.ArgumentParser("electrofit")
    sub = p.add_subparsers(dest="cmd", required=True)

    DESCRIPTIONS = {
        'step0': 'Project scaffold: copy each molecule folder from data/input/ -> process/<mol>/ (idempotent).',
        'step1': 'Initial molecule preparation: run BCC or OPT Gaussian+RESP pipeline once per molecule (produces <mol>.acpype).',
        'step2': 'Simulation setup: create run_gmx_simulation folders, stage GMX topology + MDP templates, build run.json manifest.',
        'step3': 'Production MD (GROMACS): build box, solvate, ions, EM/NVT/NPT, production trajectory per molecule.',
        'step4': 'Conformer extraction: sample frames from trajectories into extracted_conforms/ (linear|random|maxmin strategies).',
        'step5': 'Conformer charge pipeline: Gaussian/RESP per conformer (parallelisable) producing ensemble RESP charges.',
        'step6': 'Charge aggregation: compute ensemble averages, optional symmetry/group averaging, update MOL2 & acpype (experimental --remove-outlier).',
        'step7': 'Placeholder / downstream analysis stage (reserved – implementation pending).',
        'step8': 'Placeholder / export or reporting stage (reserved – implementation pending).',
    }
    def add_cmd(name: str):
        sp = sub.add_parser(name, help=DESCRIPTIONS.get(name,''), description=DESCRIPTIONS.get(name,''))
        sp.add_argument("--project", required=True, help="Path to project/case folder")
        sp.add_argument("--config", help="Optional path to electrofit.toml")
        return sp

    for cmd in STEP_MODULES:
        add_cmd(cmd)

    args, rest = p.parse_known_args()
    # Flag-only validation: if user mistakenly supplies a value after --isolate-conformer
    if args.cmd == 'step5' and '--isolate-conformer' in sys.argv:
        # finde Position
        for i, tok in enumerate(sys.argv):
            if tok == '--isolate-conformer':
                # if next token does not start with '-' and is not a subcommand -> warning
                if i+1 < len(sys.argv) and not sys.argv[i+1].startswith('-'):
                    print("[step5][warn] '--isolate-conformer' is a flag-only option; ignored value: '"+sys.argv[i+1]+"'", file=sys.stderr)
                break

    # Ensure project scaffold (create <project>/process if missing)
    project_root, _process_dir = _ensure_project_scaffold(args.project)

    # Resolve config path: explicit only (no legacy env fallback anymore).
    if args.config:
        config_path = Path(args.config).resolve()
    else:
        config_path = None

    logging.info(f"Dispatching {args.cmd} | project={project_root}"
                 + (f" | config={config_path}" if config_path else " | config=<none>"))

    _run_module(STEP_MODULES[args.cmd], project_root, config_path, rest)

if __name__ == "__main__":
    main()