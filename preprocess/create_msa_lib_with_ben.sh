#!/bin/bash
#SBATCH --job-name=filter_af3_dbs
#SBATCH --time=6:00:00           # Give it plenty of time
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8        # Script uses multiprocessing
#SBATCH --output=filter_dbs_%j.log
#SBATCH --output=/cs/labs/dina/noabirman/slurms_outs/msa_lib/%x_%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=noa.birman@mail.huji.ac.il
#SBATCH --exclude=sm-01,sm-16,sm-02,sm-03,sm-04,sm-08

# Usage: sbatch create_msa_lib_with_ben.sh <database_dir> <jsons_folder> <output_folder>

if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch $0 <database_dir> <jsons_folder> <output_folder>"
    echo ""
    echo "Arguments:"
    echo "  jsons_folder   : Path to folder containing *_data.json files"
    echo "  output_folder  : Path where filtered databases will be created"
    exit 1
fi

JSONS_FOLDER=$1
OUTPUT_FOLDER=$2

echo "Starting AlphaFold3 database filtering..."
echo "JSONs folder: $JSONS_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo ""

python3 20251215_mini_msa_db.py "$JSONS_FOLDER" "$OUTPUT_FOLDER"


echo ""
echo "Job completed at $(date)"