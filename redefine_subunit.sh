#!/bin/bash
#SBATCH --job-name=redefine_subunit

#SBATCH --output=/cs/labs/dina/noabirman/slurms_outs/redefine_subunit/redefine_subunit_%j.log
#SBATCH --mail-type=END
#SBATCH --mail-user=noa.birman@mail.huji.ac.il
#SBATCH --time=00:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

# Arguments passed when submitting the job
AF_PAIRS_DIR=$1       # required: folder with AF pairs
OPTION=${2:-}         # optional: e.g., 'ipsae'
PAE=${3:-}            # optional pae_cutoff
DIST=${4:-}           # optional dist_cutoff

# Build the python command
CMD="python3 create_graph_from_complex.py $AF_PAIRS_DIR"

if [ ! -d "$AF_PAIRS_DIR" ]; then
    echo "Error: folder $AF_PAIRS_DIR does not exist!"
    exit 1
fi

# Add optional mode (ipsae) if provided
if [ -n "$OPTION" ]; then
    CMD="$CMD $OPTION"
fi

# Add pae_cutoff if provided
if [ -n "$PAE" ]; then
    CMD="$CMD $PAE"
fi

# Add dist_cutoff if provided
if [ -n "$DIST" ]; then
    CMD="$CMD $DIST"
fi

# Print and run the command
echo "Running: $CMD"
$CMD