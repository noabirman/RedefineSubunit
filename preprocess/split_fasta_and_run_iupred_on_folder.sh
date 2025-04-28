#!/bin/bash

complex_name="$1"

# Create output directories
mkdir -p "$complex_name/input_fastas"
mkdir -p "$complex_name/iupred_outputs"

# Split multi-FASTA into one FASTA per chain
awk -v dir="$complex_name/input_fastas" '/^>/{close(out); out=dir"/"substr($0,2)".fasta"} {print > out}' "$complex_name/sequences.fasta"

# Run IUPred3 on each chain
for f in "$complex_name"/input_fastas/*.fasta; do
    fname=$(basename "$f" .fasta)
    python3 /cs/labs/dina/tsori/af3_example/RedefineSubunit/iupred3/iupred3.py "$f" long > "$complex_name/iupred_outputs/${fname}_iupred.txt"
done
  