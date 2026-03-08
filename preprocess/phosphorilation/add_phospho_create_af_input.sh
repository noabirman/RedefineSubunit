#!/bin/bash
#SBATCH --job-name=add_phospho
#SBATCH --output=/cs/labs/dina/noabirman/slurms_outs/phospho/add_phospho_%j.log
#SBATCH --mail-type=END
#SBATCH --mail-user=noa.birman@mail.huji.ac.il
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

MSA_OUTPUT_ROOT=$1
PHOSPHO_FILE=$2

echo "Running add_phospho with:"
echo "MSA output root: $MSA_OUTPUT_ROOT"
echo "Phospho file: $PHOSPHO_FILE"

python3 put_ptm_in_msa_output.py "$MSA_OUTPUT_ROOT" "$PHOSPHO_FILE"