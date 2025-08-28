#!/usr/bin/env bash
# ==============================================================================
# Electrofit IP6 paper (2025): end-to-end setup + tests + IP6 analysis runner
#
# What it does:
#   1) Clone or update the paper workspace repo
#   2) Create + activate the conda env from envs/electrofit-analysis.yml
#   3) Run integration steps (electrofit step0..step8) with your flags
#   4) Optionally run the large IP6 analysis (charges/hbonds/pp-matrix/coordination)
#
# Usage (examples):
#   ./run_electrofit.sh --ip6-project /path/to/ip6-project
#   ./run_electrofit.sh --continue-on-error
#   ./run_electrofit.sh --skip-ip6
#
# Flags:
#   --repo-url URL              : Git repo URL (default: JALK96/electrofit-ip6-paper-2025)
#   --workdir DIR               : Where to clone/work (default: directory of this script)
#   --continue-on-error         : Do not stop on the first failing step
#   --use-mamba                 : Prefer mamba for env creation if available
#   --skip-clone                : Assume repo already present; skip git clone/pull
#   --skip-env                  : Skip conda env creation/activation
#   --skip-integration          : Skip electrofit step0..step8 integration test
#   --skip-ip6                  : Skip the large IP6 analysis
#   --ip6-project DIR           : Path to the *source* IP6 project to copy into repo
#                                 (required unless --skip-ip6 is used)
#   --help                      : Print this help
#
# Notes:
#   - This script expects Conda installed and available in your shell.
#   - The conda env name is read from the YAML (fallback: "electrofit-analysis").
#   - The large IP6 copy can use substantial disk space; delete after use.
# ==============================================================================

set -Eeuo pipefail

### ---- Logging helpers --------------------------------------------------------
ts() { date +"%Y-%m-%d %H:%M:%S"; }
info() { echo "[INFO  $(ts)] $*"; }
warn() { echo "[WARN  $(ts)] $*" >&2; }
err()  { echo "[ERROR $(ts)] $*" >&2; }

### ---- Defaults ---------------------------------------------------------------
REPO_URL="https://github.com/JALK96/electrofit-ip6-paper-2025.git"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKDIR="$SCRIPT_DIR"
CONTINUE_ON_ERROR=0
USE_MAMBA=0
SKIP_CLONE=0
SKIP_ENV=0
SKIP_INTEGRATION=0
SKIP_IP6=0
IP6_PROJECT_SRC=""

### ---- Parse CLI --------------------------------------------------------------
usage() {
  sed -n '1,100p' "$0" | sed -n '1,80p' | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo-url)           REPO_URL="$2"; shift 2 ;;
    --workdir)            WORKDIR="$2"; shift 2 ;;
    --continue-on-error)  CONTINUE_ON_ERROR=1; shift ;;
    --use-mamba)          USE_MAMBA=1; shift ;;
    --skip-clone)         SKIP_CLONE=1; shift ;;
    --skip-env)           SKIP_ENV=1; shift ;;
    --skip-integration)   SKIP_INTEGRATION=1; shift ;;
    --skip-ip6)           SKIP_IP6=1; shift ;;
    --ip6-project)        IP6_PROJECT_SRC="$2"; shift 2 ;;
    --help|-h)            usage; exit 0 ;;
    *)
      err "Unknown option: $1"; usage; exit 2 ;;
  esac
done

### ---- Derived paths ----------------------------------------------------------
REPO_DIR="$(basename "$REPO_URL" .git)"                # -> electrofit-ip6-paper-2025
EF_ROOT="$WORKDIR/$REPO_DIR"                           # repo root
WS_DIR="$EF_ROOT/electrofit-workspace"                 # workspace dir inside repo
ENV_YML="$WS_DIR/envs/electrofit-analysis.yml"         # env file
INTEGRATION_PROJECT="$WS_DIR/tests/integration"        # integration tests root
INTEGRATION_CACHE="$INTEGRATION_PROJECT/cache"         # gaussian cache
# Destination inside repo where we copy the big IP6 project directory:
IP6_PROJECT_DEST="$EF_ROOT/ip6-project"

### ---- Error policy for step runner ------------------------------------------
run_step() {
  echo ">>> Running: $*"
  if "$@"; then
    echo "------------------------------------"
  else
    status=$?
    echo "!!! Step failed with exit code $status"
    echo "------------------------------------"
    if [[ $CONTINUE_ON_ERROR -eq 0 ]]; then
      exit $status
    fi
  fi
}

### ---- Conda helpers ----------------------------------------------------------
conda_hook() {
  # Make 'conda activate' work in non-interactive shell
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
  else
    err "Conda not found in PATH. Ensure conda is installed and available."
    exit 1
  fi
}

detect_env_name() {
  local yml="$1"
  local name
  if [[ -f "$yml" ]]; then
    name="$(grep -E '^\s*name\s*:\s*' "$yml" | head -n1 | sed -E 's/^\s*name\s*:\s*//')"
    if [[ -n "$name" ]]; then
      echo "$name"
      return 0
    fi
  fi
  echo "electrofit-analysis"
}

create_env_from_yml() {
  local yml="$1"
  local use_mamba="$2"
  if [[ ! -f "$yml" ]]; then
    err "Environment YAML not found: $yml"
    exit 1
  fi
  if [[ $use_mamba -eq 1 && $(command -v mamba >/dev/null; echo $?) -eq 0 ]]; then
    info "Creating env with mamba from $yml"
    mamba env create -f "$yml" || true
  else
    info "Creating env with conda from $yml"
    conda env create -f "$yml" || true
  fi
}

### ---- Git clone/update -------------------------------------------------------
clone_or_update() {
  if [[ $SKIP_CLONE -eq 1 ]]; then
    info "Skipping clone step as requested (--skip-clone)."
    return 0
  fi
  mkdir -p "$WORKDIR"
  cd "$WORKDIR"
  if [[ -d "$EF_ROOT/.git" ]]; then
    info "Repo already present: $EF_ROOT"
    cd "$EF_ROOT"
    info "Fetching latest changes..."
    git pull --rebase --autostash
  else
    info "Cloning $REPO_URL into $WORKDIR"
    git clone "$REPO_URL"
  fi
}

### ---- Integration test steps -------------------------------------------------
run_integration() {
  if [[ $SKIP_INTEGRATION -eq 1 ]]; then
    info "Skipping integration tests (--skip-integration)."
    return 0
  fi
  if [[ ! -d "$INTEGRATION_PROJECT" ]]; then
    err "Integration project directory not found: $INTEGRATION_PROJECT"
    exit 1
  fi
  info "Running integration steps in: $INTEGRATION_PROJECT"

  # Steps 0..8, with special flags where needed
  run_step electrofit step0 --project "$INTEGRATION_PROJECT"
  run_step electrofit step1 --project "$INTEGRATION_PROJECT"
  run_step electrofit step2 --project "$INTEGRATION_PROJECT"
  run_step electrofit step3 --project "$INTEGRATION_PROJECT"
  run_step electrofit step4 --project "$INTEGRATION_PROJECT" --sample 2
  run_step env ELECTROFIT_DEBUG_GAUSSIAN_CACHE="$INTEGRATION_CACHE" electrofit step5 --project "$INTEGRATION_PROJECT"
  run_step electrofit step6 --project "$INTEGRATION_PROJECT"
  run_step electrofit step7 --project "$INTEGRATION_PROJECT"
  run_step electrofit step8 --project "$INTEGRATION_PROJECT"
}

### ---- Large IP6 analysis -----------------------------------------------------
prepare_ip6_project() {
  if [[ $SKIP_IP6 -eq 1 ]]; then
    info "Skipping IP6 analysis (--skip-ip6)."
    return 1
  fi
  if [[ -z "$IP6_PROJECT_SRC" ]]; then
    err "You must provide --ip6-project /path/to/ip6-project (or use --skip-ip6)."
    exit 1
  fi
  if [[ ! -d "$IP6_PROJECT_SRC" ]]; then
    err "Given IP6 project path does not exist: $IP6_PROJECT_SRC"
    exit 1
  fi
  if [[ -d "$IP6_PROJECT_DEST" ]]; then
    warn "Destination IP6 project already exists at $IP6_PROJECT_DEST — reusing."
    return 0
  fi
  info "Copying IP6 project into repo (this may take significant disk space):"
  info "  from: $IP6_PROJECT_SRC"
  info "    to: $IP6_PROJECT_DEST"
  cp -r "$IP6_PROJECT_SRC" "$IP6_PROJECT_DEST"
}

run_ip6_analysis() {
  if [[ $SKIP_IP6 -eq 1 ]]; then
    return 0
  fi
  if [[ ! -d "$IP6_PROJECT_DEST" ]]; then
    err "IP6 project destination missing: $IP6_PROJECT_DEST"
    exit 1
  fi
  info "Running IP6 analysis on: $IP6_PROJECT_DEST"
  run_step ip6-analysis charges      --project "$IP6_PROJECT_DEST"
  run_step ip6-analysis hbonds       --project "$IP6_PROJECT_DEST"
  run_step ip6-analysis pp-matrix    --project "$IP6_PROJECT_DEST"
  run_step ip6-analysis coordination --project "$IP6_PROJECT_DEST"
}

### ---- Main -------------------------------------------------------------------
main() {
  info "Workdir     : $WORKDIR"
  info "Repo URL    : $REPO_URL"
  info "Repo root   : $EF_ROOT"

  clone_or_update
  cd "$EF_ROOT"

  if [[ $SKIP_ENV -eq 0 ]]; then
    conda_hook
    ENV_NAME="$(detect_env_name "$ENV_YML")"
    info "Conda env   : $ENV_NAME (from $ENV_YML)"
    create_env_from_yml "$ENV_YML" "$USE_MAMBA"
    info "Activating conda env: $ENV_NAME"
    conda activate "$ENV_NAME"
  else
    warn "Skipping env creation/activation (--skip-env). Make sure 'electrofit' and 'ip6-analysis' are on PATH."
  fi

  # Sanity: workspace directory should exist
  if [[ ! -d "$WS_DIR" ]]; then
    err "Workspace directory not found: $WS_DIR"
    exit 1
  fi

  run_integration
  if prepare_ip6_project; then
    run_ip6_analysis
  fi

  info "All requested tasks completed."
}

main "$@"