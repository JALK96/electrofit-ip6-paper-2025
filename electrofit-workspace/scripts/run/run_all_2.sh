#!/bin/bash

PROJECT="/home/johannal96/Publications.nobackup/2025/electrofit-ip6-paper-2025/electrofit-workspace/tests/integration"
CACHE="/home/johannal96/Publications.nobackup/2025/electrofit-ip6-paper-2025/electrofit-workspace/tests/integration/cache"
MOLECULE_NAME="IP_011101"

run_step() {
    echo ">>> Running: $*"
    "$@"
    status=$?
    if [ $status -ne 0 ]; then
        echo "!!! Step failed with exit code $status"
    fi
    echo "------------------------------------"
}

#run_step electrofit step0 --project "$PROJECT" --molecule "$MOLECULE_NAME"
#run_step electrofit step1 --project "$PROJECT" --molecule "$MOLECULE_NAME"
#run_step electrofit step2 --project "$PROJECT" --molecule "$MOLECULE_NAME"
run_step electrofit step3 --project "$PROJECT" --molecule "$MOLECULE_NAME"
run_step electrofit step4 --project "$PROJECT" --sample 2 --molecule "$MOLECULE_NAME"
run_step env ELECTROFIT_DEBUG_GAUSSIAN_CACHE="$CACHE" electrofit step5 --project "$PROJECT" --molecule "$MOLECULE_NAME"
run_step electrofit step6 --project "$PROJECT" --remove-outlier --plot-histograms --molecule "$MOLECULE_NAME"
run_step electrofit step7 --project "$PROJECT" --molecule "$MOLECULE_NAME"
run_step electrofit step8 --project "$PROJECT" --molecule "$MOLECULE_NAME"
