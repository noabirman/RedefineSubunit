#!/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --killable
#SBATCH --requeue

cd /cs/labs/dina/tsori/af3_example
source RedefineSubunit/my_venv/bin/activate
python3 CombFold/scripts/run_on_pdbs.py mll4/af_pairs/subunits_info.json mll4/pdbs mll4/comfold/
