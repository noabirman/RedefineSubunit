#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4-00:00:00

#SBATCH --mail-type=END
#SBATCH --mail-user=tsori.kislev@gmail.com

#SBATCH --exclude=sm-01,sm-16,sm-02,sm-03,sm-04,sm-08

#SBATCH --output=/cs/labs/dina/tsori/af3_example/slurms_outs/msa/%j.out
#SBATCH --error=/cs/labs/dina/tsori/af3_example/slurms_outs/msa/err/%j.out

#export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

# Documentation:
# This script runs the MSA (Multiple Sequence Alignment) process for a given input directory.
# Usage:
#   ./msa.sh <INPUT_DIR> [MAPPING_JSON] [SUBUNITS_INFO_JSON]
# Arguments:
#   <INPUT_DIR>: Path to the input directory containing JSON files.
#   [MAPPING_JSON]: (Optional) Path to the mapping JSON file. If not specified, assumes 'mapping.json' in the parent directory of INPUT_DIR.
#   [SUBUNITS_INFO_JSON]: (Optional) Path to the subunits info JSON file. If not specified, assumes 'subunits_info.json' in the parent directory of INPUT_DIR.
# Output:
#   The MSA results will be saved in a directory named 'msa_output' in the same parent directory as the input directory.
#   Pairwise MSA files will be generated using the msa_to_pairwise.py script.

export PATH="/sci/labs/dina/bshor/projects/af_combdock/tools/conda_install/miniconda3/bin:$PATH"
. "/sci/labs/dina/bshor/projects/af_combdock/tools/conda_install/miniconda3/etc/profile.d/conda.sh"
conda activate /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/alphafold3-conda

$INPUT = "$1"
# Determine parent directory and set default paths
PARENT_DIR=$(dirname "$INPUT")
$OUTPUT_DIR="${PARENT_DIR}/msa_output"

# Print the paths for debugging
echo Running MSA on : "$INPUT"
echo Output directory: "$OUTPUT_DIR"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"



python /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/alphafold3/run_alphafold.py \
  --jackhmmer_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/jackhmmer \
  --db_dir /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/databases \
  --model_dir /cs/usr/bshor/sci/installations/af3_variations/deepmind/models \
  --hmmbuild_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/hmmbuild \
  --hmmsearch_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/hmmsearch \
  --norun_inference \
  --output_dir "$OUTPUT_DIR" \
  --json_path "$INPUT"

# cd /cs/labs/dina/tsori/af3_example/RedefineSubunit/
# python preprocess/msa_to_pairwise.py "$OUTPUT_DIR" "$MAPPING_JSON" "$SUBUNITS_INFO_JSON"
