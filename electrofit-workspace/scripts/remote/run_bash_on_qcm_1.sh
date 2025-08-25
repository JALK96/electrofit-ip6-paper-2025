#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# -----------------------------
# Configuration and Setup
# -----------------------------

# Bestimme das Verzeichnis, in dem dieses Skript liegt
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Name of the Bash script to execute
BASH_SCRIPT="/home/johannal96/Publications.nobackup/2025/electrofit-ip6-paper-2025/electrofit-workspace/scripts/run/run_all_1.sh"
echo "Bash Script: $BASH_SCRIPT"

# ------------ ---------------- ------------- -----------
# Remote machine details
REMOTE_HOST="qcm04"
echo "Remote Host: $REMOTE_HOST"
# ------------ ---------------- ------------- -----------

# Generate a unique Screen session name to prevent conflicts
SCREEN_SESSION="$(basename "$(dirname "$SCRIPT_DIR")")_$(date +%Y%m%d%H%M%S)"
echo "Screen Session Name: $SCREEN_SESSION"

# Name of the conda environment to use
CONDA_ENV="electrofit"
echo "Conda Environment: $CONDA_ENV"

# -----------------------------
# Logging Setup
# -----------------------------

# Define the log file (placed in the same directory as the script)
LOG_FILE="$SCRIPT_DIR/process.log"

# Function to log messages with timestamps
log() {
    local MESSAGE="$1"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $MESSAGE" | tee -a "$LOG_FILE"
}

# Initialize the log file
touch "$LOG_FILE"
log "=== Script Started ==="

# Log initial configurations
log "Script Directory: $SCRIPT_DIR"
log "Bash Script Path: $BASH_SCRIPT"
log "Remote Host: $REMOTE_HOST"
log "Screen Session Name: $SCREEN_SESSION"

# -----------------------------
# Function Definitions
# -----------------------------

# Function to execute SSH command
execute_ssh_command() {
    local HOST="$1"
    local CMD="$2"

    log "Connecting to $HOST to execute the command."
    ssh "$HOST" "$CMD"
    SSH_STATUS=$?
    
    if [ $SSH_STATUS -ne 0 ]; then
        log "Error: SSH command failed with status $SSH_STATUS."
        exit 1
    else
        log "SSH command executed successfully on $HOST."
    fi
}

# -----------------------------
# Build the SSH Command
# -----------------------------

SSH_COMMAND="cd \"$SCRIPT_DIR\" && \
screen -dmS \"$SCREEN_SESSION\" bash -lc '
  source \"$HOME/miniconda3/etc/profile.d/conda.sh\";
  conda activate $CONDA_ENV;
  echo \"[INFO] Using python: \$(which python)\";
  echo \"[INFO] electrofit path test:\";
  python -c \"import electrofit, sys; print(\\\"electrofit from:\\\", electrofit.__file__, \\\"python:\\\", sys.executable)\";
  bash \"$BASH_SCRIPT\";
  exec bash
'"

log "Built SSH Command: $SSH_COMMAND"

# -----------------------------
# Execute the SSH Command
# -----------------------------

execute_ssh_command "$REMOTE_HOST" "$SSH_COMMAND"

# -----------------------------
# Post-Execution Confirmation
# -----------------------------

log "Successfully started '$BASH_SCRIPT' in Screen session '$SCREEN_SESSION' on $REMOTE_HOST in conda environment $CONDA_ENV."

# -----------------------------
# End of the Script
# -----------------------------

log "=== Script Completed Successfully ==="