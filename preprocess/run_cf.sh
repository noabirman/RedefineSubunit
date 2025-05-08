#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4-00:00:00

#SBATCH --mail-type=END
#SBATCH --mail-user=tsori.kislev@gmail.com

#SBATCH --output=/cs/labs/dina/tsori/af3_example/slurms_outs/cf/%j.out


# Check for required arguments
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <COMPLEX_DIR> [MODELS_DIR] [SUBUNITS_INFO_JSON] "
  exit 1
fi





cd /cs/labs/dina/tsori/af3_example
source RedefineSubunit/my_venv/bin/activate
python3 CombFold/scripts/run_on_pdbs.py mll4/af_pairs/subunits_info.json mll4/pdbs mll4/comfold/
