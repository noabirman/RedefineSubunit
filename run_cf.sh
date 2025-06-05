#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4-00:00:00

#SBATCH --mail-type=END
#SBATCH --mail-user=tsori.kislev@gmail.com

#SBATCH --output=/cs/labs/dina/tsori/af3_example/slurms_outs/CF/%j.out

set -e  # Exit on error

# Usage: ./run_cf.sh <COMPLEX_DIR> <MODE>
# MODE:
#   1: Default
#   2: High (uses high_subunits_info.json)
#   3: Trivial (uses subunits_info or iupred + trivial dirs)
#   4: Us-Trivial (subunits from combfold, trivial models/results)
#   5: All (subunits from combfold, all models, results to combfold_all)

if [ "$#" -lt 2 ]; then
  printf 'Usage: %s <COMPLEX_DIR> <MODE (1–5)>\n\
    1: Default\n\
    2: High (uses high_subunits_info.json)\n\
    3: Trivial (uses subunits_info or iupred + trivial dirs)\n\
    4: Us-Trivial (subunits from combfold, trivial models/results)\n\
    5: All (subunits from combfold,all models, results to combfold_all)\n' "$0"

  exit 1
fi

COMPLEX_DIR="$1"
MODE="$2"

case "$MODE" in
  1)
    SUBUNITS_INFO_JSON="$COMPLEX_DIR/combfold/subunits_info.json"
    MODELS_DIR="$COMPLEX_DIR/combfold/models"
    RESULTS_DIR="$COMPLEX_DIR/combfold/results"
    ;;
  2)
    SUBUNITS_INFO_JSON="$COMPLEX_DIR/combfold/high_subunits_info.json"
    MODELS_DIR="$COMPLEX_DIR/combfold/models"
    RESULTS_DIR="$COMPLEX_DIR/combfold/results_high"
    ;;
  3)
    if [ -f "$COMPLEX_DIR/iupred_subunits_info.json" ]; then
      SUBUNITS_INFO_JSON="$COMPLEX_DIR/iupred_subunits_info.json"
    else
      SUBUNITS_INFO_JSON="$COMPLEX_DIR/subunits_info.json"
    fi
    MODELS_DIR="$COMPLEX_DIR/combfold_trivial/models"
    RESULTS_DIR="$COMPLEX_DIR/combfold_trivial/results"
    ;;
  4)
    SUBUNITS_INFO_JSON="$COMPLEX_DIR/combfold/subunits_info.json"
    MODELS_DIR="$COMPLEX_DIR/combfold_trivial/models"
    RESULTS_DIR="$COMPLEX_DIR/combfold_us_trivial/results"
    ;;
  5)
    SUBUNITS_INFO_JSON="$COMPLEX_DIR/combfold/subunits_info.json"
    MODELS_DIR="$COMPLEX_DIR/combfold_all/models"
    RESULTS_DIR="$COMPLEX_DIR/combfold_all/results"
    ;;
  *)
    echo "Invalid mode: $MODE. Choose 1–5."
    exit 1
    ;;
esac

# Check inputs
if [ ! -d "$COMPLEX_DIR" ]; then
  echo "Error: Complex directory '$COMPLEX_DIR' does not exist."
  exit 1
fi

if [ ! -f "$SUBUNITS_INFO_JSON" ]; then
  echo "Error: Subunits info JSON file '$SUBUNITS_INFO_JSON' does not exist."
  exit 1
fi

if [ ! -d "$MODELS_DIR" ]; then
  echo "Error: Models directory '$MODELS_DIR' does not exist."
  exit 1
fi

mkdir -p "$RESULTS_DIR"

cd /cs/labs/dina/tsori/af3_example
source RedefineSubunit/my_venv/bin/activate

echo "Running CombFold mode $MODE"
echo "Subunits:   $SUBUNITS_INFO_JSON"
echo "Models dir: $MODELS_DIR"
echo "Results to: $RESULTS_DIR"
echo -------------------------------------------------------------------------------------

python3 CombFold/scripts/run_on_pdbs.py "$SUBUNITS_INFO_JSON" "$MODELS_DIR" "$RESULTS_DIR"

echo
echo "Done. Output log: $RESULTS_DIR/_unified_representation/assembly_output/output.log"
