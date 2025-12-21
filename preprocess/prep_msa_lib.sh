#!/bin/bash
#SBATCH --job-name=af3_filter
#SBATCH --output=%x_%j.out    # Output logs
#SBATCH --output=/cs/labs/dina/noabirman/slurms_outs/msa_lib/%x_%j.out
#SBATCH --mem=16G             # Memory per job
#SBATCH --cpus-per-task=5     # CPU cores
#SBATCH --time=4:00:00        # Max time limit (2 hours per folder)

#SBATCH --mail-type=END
#SBATCH --mail-user=noa.birman@mail.huji.ac.il
#SBATCH --exclude=sm-01,sm-16,sm-02,sm-03,sm-04,sm-08

# 1. Get the folder name from the command line argument
TARGET_DIR=$1

# 2. Go into the folder
cd "$TARGET_DIR" || exit

# 3. Run the python script (using full path)
python3 /cs/labs/dina/noabirman/RedefineSubunit/preprocess/extract_MSA_templates.py