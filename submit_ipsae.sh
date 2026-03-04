#!/bin/bash
#SBATCH --job-name=ipsae_score

#SBATCH --output=/cs/labs/dina/noabirman/slurms_outs/ipsae/ipsae_%j.log
#SBATCH --mail-type=END
#SBATCH --mail-user=noa.birman@mail.huji.ac.il
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

ROOT_DIR=$1
PAE=$2
DIST=$3

echo "Running IPSAE with:"
echo "Root dir: $ROOT_DIR"
echo "PAE cutoff: $PAE"
echo "Distance cutoff: $DIST"

python3 run_ipsae_on_af_pairs.py "$ROOT_DIR" "$PAE" "$DIST"