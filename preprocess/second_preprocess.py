#!/usr/bin/env python3
# changed to create pairwise msa input instead if single
"""
Reads a subunits_info JSON, assigns alphabetical labels, and creates 
pairwise AlphaFold 3 JSONs. Routes pairs to msa_inputs or msa_output 
based on sequence lengths.
"""
import os
import json
import string
import itertools
import argparse

def combine_and_pair_subunits(input_file, dir_path, min_length):
    # Load the input JSON file
    with open(input_file, 'r') as f:
        subunits_data = json.load(f)

    # Create directories
    msa_inputs_dir = os.path.join(dir_path, "msa_inputs")
    msa_outputs_dir = os.path.join(dir_path, "msa_output")
    os.makedirs(msa_inputs_dir, exist_ok=True)
    os.makedirs(msa_outputs_dir, exist_ok=True)

    mapping = {}
    chains = {}

    # Generate labels A, B, ..., Z, AA, AB...
    def generate_labels():
        alphabet = string.ascii_uppercase
        for letter in alphabet:
            yield letter
        for pair in itertools.product(alphabet, repeat=2):
            yield ''.join(pair)

    label_generator = generate_labels()

    # Step 1: Parse subunits and determine length status
    for subunit_name, subunit_info in subunits_data.items():
        new_label = next(label_generator)
        seq = subunit_info["sequence"]
        seq_len = len(seq)
        start_res = subunit_info["start_res"]
        
        # Save mapping info
        mapping[new_label] = {
            "chain_id": subunit_info["chain_names"][0],
            "start": start_res,
            "end": start_res + seq_len - 1
        }
        
        # Store sequence and short/long status for combinatorial pairing
        chains[new_label] = {
            "sequence": seq,
            "is_short": seq_len < min_length
        }

    # Step 2: Generate all pairwise combinations
    chain_ids = sorted(chains.keys())
    print(f"Chains detected: {chain_ids}\n")

    for A, B in itertools.combinations(chain_ids, 2):
        pair_name = f"{A}_{B}"
        is_A_short = chains[A]["is_short"]
        is_B_short = chains[B]["is_short"]
        
        # Base JSON structure for an AF3 pair
        pair_data = {
            "name": pair_name,
            "modelSeeds": [1],
            "sequences": [
                {"protein": {"id": A, "sequence": chains[A]["sequence"]}},
                {"protein": {"id": B, "sequence": chains[B]["sequence"]}}
            ],
            "dialect": "alphafold3",
            "version": 1
        }

        # Step 3: Route based on length logic
        if is_A_short and is_B_short:
            # Short + Short: Bypass MSA pipeline completely
            pair_data["sequences"][0]["protein"].update({"unpairedMsa": "", "pairedMsa": "", "templates": []})
            pair_data["sequences"][1]["protein"].update({"unpairedMsa": "", "pairedMsa": "", "templates": []})
            
            output_subdir = os.path.join(msa_outputs_dir, pair_name.lower())
            os.makedirs(output_subdir, exist_ok=True)
            out_path = os.path.join(output_subdir, f"{pair_name.lower()}_data.json")
            print(f"Short+Short Pair: Saved to msa_output -> {out_path}")
        else:
            # Long+Long OR Long+Short: Sent to inputs for MSA processing
            out_path = os.path.join(msa_inputs_dir, f"{pair_name}.json")
            print(f"Needs MSA Pipeline: Saved to msa_inputs -> {out_path}")

        with open(out_path, 'w') as f:
            json.dump(pair_data, f, indent=4)

    # Save mapping file
    mapping_file_path = os.path.join(dir_path, "chain_id_mapping.json")
    with open(mapping_file_path, 'w') as f:
        json.dump(mapping, f, indent=4)
    print(f"\nMapping file saved to {mapping_file_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process subunits and create pairwise AF3 JSONs.")
    parser.add_argument("combfold_dir", help="Directory path for outputs")
    parser.add_argument("--min_length", type=int, default=10, help="Minimum sequence length for MSA (default: 10)")
    parser.add_argument("--input_file", help="Path to subunits_info.json (optional)")
    
    args = parser.parse_args()
    
    # Default to subunits_info.json in the combfold_dir if not explicitly provided
    input_file = args.input_file if args.input_file else os.path.join(os.path.abspath(args.combfold_dir), 'subunits_info.json')
    
    print(f"Directory path: {os.path.abspath(args.combfold_dir)}")
    print(f"Input file: {input_file}")
    print(f"Minimum sequence length for MSA: {args.min_length}\n")
    
    combine_and_pair_subunits(input_file, os.path.abspath(args.combfold_dir), args.min_length)