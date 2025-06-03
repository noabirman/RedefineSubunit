#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4-00:00:00

#SBATCH --mail-type=END
#SBATCH --mail-user=tsori.kislev@gmail.com

#SBATCH --output=/cs/labs/dina/tsori/af3_example/slurms_outs/CF/%j.out

set -e  # Exit immediately if any command fails

# Documentation:
# This script runs the CombFold pipeline on a specified complex directory.
# Usage:
#   ./run_cf.sh --complex_dir <COMPLEX_DIR> [--high_only] [--subunits <SUBUNITS_INFO_JSON>] [--models <MODELS_DIR>]
# Arguments:
#   --complex_dir <COMPLEX_DIR>         Path to the complex directory (required)
#   --high_only                        (Optional) If present, only use subunits with '_high_' in their key.
#   --subunits <SUBUNITS_INFO_JSON>     (Optional) Path to the subunits info JSON file.
#   --models <MODELS_DIR>               (Optional) Path to the models directory.
# Output:
#   Results will be saved in <COMPLEX_DIR>/combfold/results.

# Default values
COMPLEX_DIR=""
HIGH_ONLY=0
SUBUNITS_INFO_JSON=""
MODELS_DIR=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --complex_dir)
      COMPLEX_DIR="$2"
      shift 2
      ;;
    --high_only)
      HIGH_ONLY=1
      shift
      ;;
    --subunits)
      SUBUNITS_INFO_JSON="$2"
      shift 2
      ;;
    --models)
      MODELS_DIR="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      echo "Usage: $0 --complex_dir <COMPLEX_DIR> [--high_only] [--subunits <SUBUNITS_INFO_JSON>] [--models <MODELS_DIR>]"
      exit 1
      ;;
  esac
done

# Check required argument
if [ -z "$COMPLEX_DIR" ]; then
  echo "Error: --complex_dir is required."
  exit 1
fi

# Set defaults if not provided
if [ -z "$SUBUNITS_INFO_JSON" ]; then
  SUBUNITS_INFO_JSON="$COMPLEX_DIR/combfold/subunits_info.json"
fi
if [ -z "$MODELS_DIR" ]; then
  MODELS_DIR="$COMPLEX_DIR/combfold/models"
fi

# Validate COMPLEX_DIR
if [ ! -d "$COMPLEX_DIR" ]; then
  echo "Error: Complex directory '$COMPLEX_DIR' does not exist."
  exit 1
fi

# Validate SUBUNITS_INFO_JSON
if [ ! -f "$SUBUNITS_INFO_JSON" ]; then
  echo "Error: Subunits info JSON file '$SUBUNITS_INFO_JSON' does not exist."
  exit 1
fi

# If high_only is set, filter the subunits_info.json
if [ "$HIGH_ONLY" -eq 1 ]; then
  python3 RedefineSubunit/filter_high_subunits.py "$SUBUNITS_INFO_JSON"
  SUBUNITS_INFO_JSON="$(dirname "$SUBUNITS_INFO_JSON")/high_subunits_info.json"
fi

RESULTS_DIR="$COMPLEX_DIR/combfold/results"
rm -rf "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

cd /cs/labs/dina/tsori/af3_example
source RedefineSubunit/my_venv/bin/activate

echo "running cif_to_pdb.py.."
echo
echo -------------------------------------------------------------------------------------
echo

python3 RedefineSubunit/cif_to_pdb.py "$COMPLEX_DIR"

echo "Running CombFold on $COMPLEX_DIR:"
echo
echo -------------------------------------------------------------------------------------
echo
python3 CombFold/scripts/run_on_pdbs.py "$SUBUNITS_INFO_JSON" "$MODELS_DIR" "$RESULTS_DIR"

echo
echo  "output log path: $RESULTS_DIR/_unified_representation/assembly_output/output.log"