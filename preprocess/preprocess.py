"""
Script Name: preprocess.py

Description:
    Preprocessing pipeline for AlphaFold 3 input:
    1) Converts AF3 input JSON to a multi-FASTA file.
    2) Splits the multi-FASTA into individual FASTA files (one per chain).
    3) Runs IUPred3 on each sequence.
    4) Splits long sequences (>1000 residues) only at disordered regions (avoiding ordered blocks).
    5) Renames chain IDs to alphabetical characters and saves a mapping.
    6) Creates a cleaned JSON file ready for AlphaFold 3 MSA input.
    7) splits the JSON file into individual JSON files (one per chain).

Usage:
    python preprocess.py <input_json> <complex_name>

Author: Noa Birman
Date: 2025-04-21
"""

import sys
import subprocess
import os
import json
from rename_protein_ids import rename_protein_ids

def json_to_fasta(json_file:str, complexname:str): #should be path/synapse.json
    # Load your JSON file
    with open(json_file) as f:
        data = json.load(f)

    # Write to a FASTA file
    with open(f"{complexname}/sequences.fasta", "w") as fasta:
        for entry in data["sequences"]:
            pid = entry["protein"]["id"]
            seq = entry["protein"]["sequence"]
            fasta.write(f">{pid}\n{seq}\n")

def subunit_to_fasta(json_file: str, complexname: str):
    """
    Converts a subunits-style JSON file to a FASTA file.

    Args:
        json_file (str): Path to the subunits JSON file.
        complexname (str): Directory to save the FASTA file.
    """
    # Load the JSON file
    with open(json_file) as f:
        data = json.load(f)

    # Copy the input JSON file to the output directory with the name 'subunits.json'
    subunits_json_path = os.path.join(complexname, "subunits_info.json")
    with open(subunits_json_path, "w") as f:
        json.dump(data, f, indent=4)

    # Write to a FASTA file
    with open(f"{complexname}/sequences.fasta", "w") as fasta:
        for subunit in data.values():
            sequence_info = {
                "id": subunit["chain_names"][0],
                "sequence": subunit["sequence"]
            }
            fasta.write(f">{sequence_info['id']}\n{sequence_info['sequence']}\n")

def split_af3_json_by_chain(input_json_path, output_dir="msa_input"):
    """
    Splits an AF3-style JSON file into separate files per protein sequence.

    Args:
        input_json_path (str): Path to the full AF3 input JSON file.
        output_dir (str): Directory to save individual chain JSONs.
    """
    os.makedirs(output_dir, exist_ok=True)

    with open(input_json_path, 'r') as f:
        af3_data = json.load(f)

    for protein_entry in af3_data["sequences"]:
        chain_id = protein_entry["protein"]["id"]
        chain_data = {
            "name": f"{chain_id}",
            "modelSeeds": af3_data.get("modelSeeds", [1]),
            "sequences": [protein_entry],
            "dialect": af3_data.get("dialect", "alphafold3"),
            "version": af3_data.get("version", 1)
        }

        output_path = os.path.join(output_dir, f"{chain_id}.json")
        with open(output_path, 'w') as out_file:
            json.dump(chain_data, out_file, indent=2)

        print(f"Saved {output_path}")

import json
from collections import defaultdict

def merge_duplicate_chain_sequences(json_path):
    # Load subunits JSON
    with open(json_path, "r") as f:
        subunits = json.load(f)

    # Group subunits by chain name
    chain_to_subunits = defaultdict(list)
    for subunit_name, info in subunits.items():
        for chain in info["chain_names"]:
            chain_to_subunits[chain].append((subunit_name, info))

    merged_subunits = {}

    # Merge by chain, keep the name of the first subunit
    for chain_name, entries in chain_to_subunits.items():
        sorted_entries = sorted(entries, key=lambda x: x[1]["start_res"])
        first_name, first_info = sorted_entries[0]
        full_sequence = "".join(info["sequence"] for _, info in sorted_entries)

        merged_subunits[first_name] = {
            "name": first_info["name"],
            "chain_names": [chain_name],
            "start_res": first_info["start_res"],
            "sequence": full_sequence
        }

    # Save merged subunits back to the same JSON file
    with open(json_path, "w") as f:
        json.dump(merged_subunits, f, indent=4)

    print(f"Merged subunits saved to {json_path}")

if __name__ == '__main__':
    # === Input Arguments ===
    if len(sys.argv) < 2 or len(sys.argv) > 4:
        print("Usage: python3 preprocess.py <dir_path> [json_path] [mode]")
        sys.exit(1)

    # Parse arguments
    dir_path = os.path.abspath(sys.argv[1])
    json_path = os.path.join(dir_path, 'subunits_info.json') if len(sys.argv) == 2 else os.path.abspath(sys.argv[2])
    mode = sys.argv[3] if len(sys.argv) == 4 else None

    # Validate paths
    if not os.path.isdir(dir_path):
        print(f"Error: '{dir_path}' is not a valid directory.")
        sys.exit(1)

    if not os.path.isfile(json_path):
        print(f"Error: '{json_path}' does not exist or is not a valid file.")
        sys.exit(1)

    # Print the parsed arguments
    print(f"Directory Path: {dir_path}")
    print(f"JSON Path: {json_path}")
    if mode:
        print(f"Mode: {mode}")

    # === Paths ===
    script_dir = os.path.dirname(os.path.abspath(__file__))
    split_script = os.path.join(script_dir, "split_fasta_and_run_iupred_on_folder.sh")
    os.makedirs(dir_path, exist_ok=True)
    split_fasta_dir = os.path.join(dir_path, "input_fastas") # Where per-chain FASTA files go
    split_mapping_file = os.path.join(dir_path, "iupred_split_mapping.json")
    split_fasta_out = os.path.join(dir_path, "iupred_split_sequences.fasta")
    iupred_outputs_path = os.path.join(dir_path, "iupred_outputs")
    fasta_path = os.path.join(dir_path, "sequences.fasta")
    af3_json = os.path.join(dir_path, "af3_input.json")
    af3_json_renamed = os.path.join(dir_path, "af3_input_renamed.json")
    mapping_file_path = os.path.join(dir_path, "chain_id_mapping.json")
    msa_inputs_path = os.path.join(dir_path, "msa_inputs")

    # Print preprocessing message
    print(f"Preprocessing \"{dir_path}\"...")
    # === Step 1: Convert JSON to full FASTA ===
    print("Converting JSON to FASTA...")
    if mode:
        json_to_fasta(json_path, dir_path)
    else:
        merge_duplicate_chain_sequences(json_path)
        subunit_to_fasta(json_path, dir_path)
    # === Step 2: Split FASTA and run IUPred3 on each chain ===
    print("Splitting FASTA and running IUPred3...")
    subprocess.run(["bash", split_script, dir_path], check=True)

    # === Step 3: Split long sequences based on IUPred3 disorder ===
    print("Splitting long sequences at disordered regions...")
    subprocess.run([
        "python3",os.path.join(script_dir, "split_sequences_on_disorder.py"),
        iupred_outputs_path,
        fasta_path,
        dir_path
    ], check=True)

    # === Step 4: Rename chain IDs and save mapping ===
    print("Renaming chain IDs...")

    # Rename chain IDs and save mapping
    rename_protein_ids(af3_json, af3_json_renamed, mapping_file_path, split_mapping_file)

    # === Step 5: Split the renamed json into per chain===
    split_af3_json_by_chain(af3_json_renamed, msa_inputs_path)

    print("Preprocessing complete.")
