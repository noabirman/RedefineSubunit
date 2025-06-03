#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4-00:00:00

#SBATCH --mail-type=END
#SBATCH --mail-user=tsori.kislev@gmail.com

#SBATCH --output=/cs/labs/dina/tsori/af3_example/slurms_outs/CF/%j.out

set -e  # Exit immediately if any command fails

# Documentation:
# Usage:
#   ./run_cf.sh <COMPLEX_DIR> [MODE]
# Arguments:
#   <COMPLEX_DIR>: Path to the complex directory (required)
#   [MODE]: (Optional) 'high' or 'trivial'
#     - If not specified or any other value: default (combfold layout)
#     - If 'high': use filtered subunits_info.json and combfold layout
#     - If 'trivial': use trivial layout (subunits_info.json and models in COMPLEX_DIR, results in COMPLEX_DIR/results)

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <COMPLEX_DIR> [MODE]"
  exit 1
fi

COMPLEX_DIR="$1"
MODE="${2:-default}"

if [ "$MODE" = "trivial" ]; then
  SUBUNITS_INFO_JSON="$COMPLEX_DIR/subunits_info.json"
  MODELS_DIR="$COMPLEX_DIR/models"
  RESULTS_DIR="$COMPLEX_DIR/results"
  CIF_INPUT="$COMPLEX_DIR"
elif [ "$MODE" = "high" ]; then
  SUBUNITS_INFO_JSON="$COMPLEX_DIR/combfold/subunits_info.json"
  MODELS_DIR="$COMPLEX_DIR/combfold/models"
  RESULTS_DIR="$COMPLEX_DIR/combfold/results"
  CIF_INPUT="$COMPLEX_DIR/combfold"
else
  SUBUNITS_INFO_JSON="$COMPLEX_DIR/combfold/subunits_info.json"
  MODELS_DIR="$COMPLEX_DIR/combfold/models"
  RESULTS_DIR="$COMPLEX_DIR/combfold/results"
  CIF_INPUT="$COMPLEX_DIR/combfold"
fi

# Validate input directories and files
if [ ! -d "$COMPLEX_DIR" ]; then
  echo "Error: Complex directory '$COMPLEX_DIR' does not exist."
  exit 1
fi
if [ ! -d "$MODELS_DIR" ]; then
  echo "Error: Models directory '$MODELS_DIR' does not exist."
  exit 1
fi
if [ ! -f "$SUBUNITS_INFO_JSON" ]; then
  echo "Error: Subunits info JSON file '$SUBUNITS_INFO_JSON' does not exist."
  exit 1
fi

# If mode is 'high', filter subunits_info.json
if [ "$MODE" = "high" ]; then
  python3 RedefineSubunit/filter_high_subunits.py "$SUBUNITS_INFO_JSON"
  SUBUNITS_INFO_JSON="$(dirname "$SUBUNITS_INFO_JSON")/high_subunits_info.json"
fi

rm -rf "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

cd /cs/labs/dina/tsori/af3_example
source RedefineSubunit/my_venv/bin/activate

echo "running cif_to_pdb.py on $CIF_INPUT"
echo
echo -------------------------------------------------------------------------------------
echo

python3 RedefineSubunit/cif_to_pdb.py "$CIF_INPUT"

echo "Running CombFold on $COMPLEX_DIR:"
echo
echo -------------------------------------------------------------------------------------
echo
python3 CombFold/scripts/run_on_pdbs.py "$SUBUNITS_INFO_JSON" "$MODELS_DIR" "$RESULTS_DIR"

echo
echo  "output log path: $RESULTS_DIR/_unified_representation/assembly_output/output.log"