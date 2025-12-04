#!/bin/bash

#SBATCH --mem=40G
#SBATCH --time=72:00:00
#SBATCH --gres=gg:g4:1
#SBATCH --exclude=creek-01,creek-02,firth-02,firth-01

#SBATCH --mail-type=END
#SBATCH --mail-user=noa.birman@mail.huji.ac.il

#SBATCH --output=/cs/labs/dina/noabirman/slurms_outs/af/%j.out

export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

export PATH="/sci/labs/dina/bshor/projects/af_combdock/tools/conda_install/miniconda3/bin:$PATH"
. "/sci/labs/dina/bshor/projects/af_combdock/tools/conda_install/miniconda3/etc/profile.d/conda.sh"
conda activate /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/alphafold3-conda

INPUT_DIR="$1"

if [ ! -d "$INPUT_DIR" ]; then
  echo "Error: Input directory '$INPUT_DIR' does not exist."
  exit 1
fi

PARENT_DIR=$(dirname "$INPUT_DIR")
OUTPUT_DIR="$PARENT_DIR/af_pairs"
MISSING_DIR="$PARENT_DIR/tmp_missing_af_input"

mkdir -p "$OUTPUT_DIR"

rm -rf "$MISSING_DIR"
mkdir "$MISSING_DIR"

# Print the paths for debugging
echo Running AF on directory: "$INPUT_DIR"
echo Output directory: "$OUTPUT_DIR"
echo tmp directory "$MISSING_DIR"

count=0
for infile in "$INPUT_DIR"/*.json; do
  chain=$(basename "$infile" .json)
  chain=${chain,,}

  # Check if output directory exists and has any .json
  if [[ -d "$OUTPUT_DIR/$chain" ]] && compgen -G "$OUTPUT_DIR/$chain/*.json" > /dev/null; then
    continue
  fi

  cp "$infile" "$MISSING_DIR/"
  ((count++))
done

if [[ $count -eq 0 ]]; then
  echo "✓ All inputs are already processed. Exiting."
  exit 0
fi

echo "⧗ Found $count missing chains. Running AlphaFold3..."

python /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/alphafold3/run_alphafold.py \
  --jackhmmer_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/jackhmmer \
  --db_dir /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/databases \
  --model_dir /cs/usr/bshor/sci/installations/af3_variations/deepmind/models \
  --hmmbuild_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/hmmbuild \
  --hmmsearch_binary_path /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/hmmer/bin/hmmsearch \
  --norun_data_pipeline \
  --output_dir "$OUTPUT_DIR" \
  --input_dir "$MISSING_DIR" \
  --flash_attention_implementation xla

