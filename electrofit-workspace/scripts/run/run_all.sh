#!/bin/bash

PROJECT="/Users/arthurlaux/Pythonprogramms/electrofit-ip6-paper-2025/electrofit-workspace/tests/integration"
CACHE="/Users/arthurlaux/Pythonprogramms/electrofit-ip6-paper-2025/electrofit-workspace/tests/integration/cache"

run_step() {
    echo ">>> Running: $*"
    "$@"
    status=$?
    if [ $status -ne 0 ]; then
        echo "!!! Step failed with exit code $status"
    fi
    echo "------------------------------------"
}

run_step electrofit step0 --project "$PROJECT" 
run_step electrofit step1 --project "$PROJECT" 
run_step electrofit step2 --project "$PROJECT" 
run_step electrofit step3 --project "$PROJECT" 
run_step electrofit step4 --project "$PROJECT" --sample 2 
run_step env ELECTROFIT_DEBUG_GAUSSIAN_CACHE="$CACHE" electrofit step5 --project "$PROJECT" 
run_step electrofit step6 --project "$PROJECT" 
run_step electrofit step7 --project "$PROJECT" 
run_step electrofit step8 --project "$PROJECT" 
