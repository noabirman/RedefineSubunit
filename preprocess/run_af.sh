#!/bin/bash

#SBATCH --job-name="af3"
#SBATCH --mem=40G
#SBATCH --time=72:00:00
#SBATCH --gres gg:g4:1
#SBATCH --exclude=creek-01,creek-02,firth-02,firth-01

#SBATCH --mail-type=END
#SBATCH --mail-user=tsori.kislev@gmail.com

#SBATCH --output=/cs/labs/dina/tsori/af3_example/slurms_outs/AF/%j.out

export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"


# Documentation:
# This script runs the AlphaFold3 pipeline on a specified input directory.
# Usage:
#   ./run_af.sh <INPUT_DIR>
# Arguments:
#   <INPUT_DIR>: Path to the input directory containing the required files.
# Output:
#   The results will be saved in a directory named 'af_pairs' in the same parent directory as the input directory.

export PATH="/sci/labs/dina/bshor/projects/af_combdock/tools/conda_install/miniconda3/bin:$PATH"
. "/sci/labs/dina/bshor/projects/af_combdock/tools/conda_install/miniconda3/etc/profile.d/conda.sh"
conda activate /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/alphafold3-conda

if [ -z "$1" ]; then
  echo "Usage: $0 <INPUT_DIR>"
  exit 1
fi

INPUT_DIR="$1"

# Validate input directory
if [ ! -d "$INPUT_DIR" ]; then
  echo "Error: Input directory '$INPUT_DIR' does not exist."
  exit 1
fi

# Determine output directory
PARENT_DIR=$(dirname "$INPUT_DIR")
OUTPUT_DIR="$PARENT_DIR/af_pairs"

echo "Run AF3 on directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

python /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/alphafold3/run_alphafold.py \
  --jackhmmer_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/jackhmmer \
  --db_dir /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/databases \
  --model_dir /cs/usr/bshor/sci/installations/af3_variations/deepmind/models \
  --hmmbuild_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/hmmbuild \
  --hmmsearch_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/hmmsearch \
  --norun_data_pipeline \
  --output_dir "$OUTPUT_DIR" \
  --input_dir "$INPUT_DIR"\
  --flash_attention_implementation xla

