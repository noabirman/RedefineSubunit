#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4-00:00:00

#SBATCH --mail-type=END
#SBATCH --mail-user=noa.birman@mail.huji.ac.il

#SBATCH --exclude=sm-01,sm-16,sm-02,sm-03,sm-04,sm-08

#SBATCH --output=/cs/labs/dina/noabirman/slurms_outs/msa/%j.out

#export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

# Documentation:
# This script runs the MSA (Multiple Sequence Alignment) process for a given input directory.
# Usage:
#   ./msa.sh <INPUT_DIR> [DB_DIR]
# Example:
#   sbatch msa.sh msa_inputs
#   sbatch msa.sh msa_inputs /path/to/mini_msa_db

# Arguments:
#   <INPUT_DIR>: Path to the input directory containing JSON files.
#   <DB_DIR>: Path to the MSA DB dir, if mini MSA library was created.
# Output:
#   The MSA results will be saved in a directory named 'msa_output' in the same parent directory as the input directory.
#   Pairwise MSA files will be generated using the msa_to_pairwise.py script if needed.

export PATH="/sci/labs/dina/bshor/projects/af_combdock/tools/conda_install/miniconda3/bin:$PATH"
. "/sci/labs/dina/bshor/projects/af_combdock/tools/conda_install/miniconda3/etc/profile.d/conda.sh"
conda activate /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/alphafold3-conda

# Check for required arguments
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <INPUT_DIR> [DB_DIR]"
  exit 1
fi

INPUT_DIR="$1"

# Validate input directory
if [ ! -d "$INPUT_DIR" ]; then
  echo "Error: Input directory '$INPUT_DIR' does not exist."
  exit 1
fi

# Default database directory
DEFAULT_DB_DIR="/cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/databases"

# Optional override DB_DIR
if [ "$#" -ge 2 ]; then
  DB_DIR="$2"
else
  DB_DIR="$DEFAULT_DB_DIR"
fi

# Validate DB_DIR
if [ ! -d "$DB_DIR" ]; then
  echo "Error: DB_DIR '$DB_DIR' does not exist."
  exit 1
fi
echo "Using DB directory: $DB_DIR"

# Determine parent directory and set default paths
PARENT_DIR=$(dirname "$INPUT_DIR")


# Determine output directory
OUTPUT_DIR="$PARENT_DIR/msa_output"

# Ensure the output directory has a unique name
#if [ -d "$OUTPUT_DIR" ]; then
 # i=2
  #while [ -d "${OUTPUT_DIR}_${i}" ]; do
   # ((i++))
  #done
  #OUTPUT_DIR="${OUTPUT_DIR}_${i}"
#fi

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Print the paths for debugging
echo Running MSA on directory: "$INPUT_DIR"

echo Output directory: "$OUTPUT_DIR"

python /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/alphafold3/run_alphafold.py \
  --jackhmmer_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/jackhmmer \
  --db_dir "$DB_DIR" \
  --model_dir /cs/usr/bshor/sci/installations/af3_variations/deepmind/models \
  --hmmbuild_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/hmmbuild \
  --hmmsearch_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/hmmsearch \
  --norun_inference \
  --output_dir "$OUTPUT_DIR" \
  --input_dir "$INPUT_DIR"

