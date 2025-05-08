#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4-00:00:00

#SBATCH --mail-type=END
#SBATCH --mail-user=tsori.kislev@gmail.com

#SBATCH --output=/cs/labs/dina/tsori/af3_example/slurms_outs/cf/%j.out

# Check for required arguments
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <COMPLEX_DIR> [MODELS_DIR] [SUBUNITS_INFO_JSON]"
  exit 1
fi

# Get input arguments
COMPLEX_DIR="$1"
MODELS_DIR="${2:-$COMPLEX_DIR/combfold/models}"  # Default to COMPLEX_DIR/combfold/models if not provided
SUBUNITS_INFO_JSON="${3:-$COMPLEX_DIR/subunits_info.json}"  # Default to COMPLEX_DIR/subunits_info.json if not provided

# Validate COMPLEX_DIR
if [ ! -d "$COMPLEX_DIR" ]; then
  echo "Error: Complex directory '$COMPLEX_DIR' does not exist."
  exit 1
fi

# Validate MODELS_DIR
if [ ! -d "$MODELS_DIR" ]; then
  echo "Error: Models directory '$MODELS_DIR' does not exist."
  exit 1
fi

# Validate SUBUNITS_INFO_JSON
if [ ! -f "$SUBUNITS_INFO_JSON" ]; then
  echo "Error: Subunits info JSON file '$SUBUNITS_INFO_JSON' does not exist."
  exit 1
fi

# Activate the virtual environment
cd /cs/labs/dina/tsori/af3_example
source RedefineSubunit/my_venv/bin/activate

# Run the Python script
RESULTS_DIR="$COMPLEX_DIR/combfold/results"
python3 CombFold/scripts/run_on_pdbs.py "$SUBUNITS_INFO_JSON" "$MODELS_DIR" "$RESULTS_DIR"