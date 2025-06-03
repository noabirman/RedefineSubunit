#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4-00:00:00

#SBATCH --mail-type=END
#SBATCH --mail-user=tsori.kislev@gmail.com

#SBATCH --output=/cs/labs/dina/tsori/af3_example/slurms_outs/CF/%j.out

set -e  # Exit immediately if any command fails


# Check for required arguments
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <COMPLEX_DIR> [SUBUNITS_INFO_JSON] [MODELS_DIR]"
  exit 1
fi

# Documentation:
# This script runs the CombFold pipeline on a specified complex directory.
# Usage:
#   ./run_cf.sh <COMPLEX_DIR> [MODELS_DIR] [SUBUNITS_INFO_JSON]
# Arguments:
#   <COMPLEX_DIR>: Path to the complex directory containing input files.
#   [MODELS_DIR]: (Optional) Path to the models directory. Defaults to <COMPLEX_DIR>/combfold/models.
#   [SUBUNITS_INFO_JSON]: (Optional) Path to the subunits info JSON file. Defaults to <COMPLEX_DIR>/subunits_info.json.
# Output:
#   Results will be saved in <COMPLEX_DIR>/combfold/results.

# Get input arguments
COMPLEX_DIR="$1"
SUBUNITS_INFO_JSON="${2:-$COMPLEX_DIR/combfold/subunits_info.json}"  # Default to COMPLEX_DIR/subunits_info.json if not provided
MODELS_DIR="${3:-$COMPLEX_DIR/combfold/models}"  # Default to COMPLEX_DIR/combfold/ if not provided

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
# Ensure the results directory exists
RESULTS_DIR="$COMPLEX_DIR/combfold/results"

rm -rf "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

# Activate the virtual environment
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

