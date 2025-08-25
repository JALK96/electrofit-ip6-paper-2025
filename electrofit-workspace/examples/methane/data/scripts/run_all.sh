#!/bin/bash

PROJECT="electrofit-workspace/examples/methane"

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
run_step electrofit step4 --project "$PROJECT" --sample 50 
run_step electrofit step5 --project "$PROJECT" 
run_step electrofit step6 --project "$PROJECT" --plot-histograms 
run_step electrofit step7 --project "$PROJECT" 
run_step electrofit step8 --project "$PROJECT" 
