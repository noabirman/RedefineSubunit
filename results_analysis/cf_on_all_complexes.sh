#!/bin/bash

# List of complexes
complexes=(
    "7arc" "7e8t" "7oba" "7ozn" "7pkn" "7qve" "7t3b" "7uic" "7wff"
    "7xho" "7zkq" "8a3t" "8a5o" "8adl" "8adn" "8bel" "8cte" "8f5o"
    "8fee" "8hil" "7p2y" "7qru" "7t2r" "7use" "8e9g"
)


BASE_DIR = "/cs/labs/dina/noabirman/complexes/DONE_MSA2"
RUN_SCRIPT = "/cs/labs/dina/noabirman/RedefineSubunit/combfold/run_cf.sh"

for complex_name in "${complexes[@]}"; do
COMPLEX_PATH = "$BASE_DIR/$complex_name"

if [ ! -d "$COMPLEX_PATH"]; then
echo
"Complex folder not found: $COMPLEX_PATH, skipping."
continue
fi

echo
"Submitting SLURM job for $complex_name..."

# Submit the job
sbatch - -requeue - -job - name = "CF_${complex_name}" "$RUN_SCRIPT" "$COMPLEX_PATH"
1
done