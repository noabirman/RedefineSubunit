import pandas as pd
from itertools import combinations
import json
import os
#from txt_to_msa2 import convert_fasta_to_msa_json

"""
    Extract only seqs that seen in pair at "PSD_residue_pairs.csv" and creates an almost msa input.
"""
def extract_unique_pairs(csv_file, col1_name, col2_name, output_file):
    # Load the CSV file
    df = pd.read_csv(csv_file)

    # Extract rows where the two columns have different values
    df_filtered = df[df[col1_name] != df[col2_name]]

    # Create unique pairs (sorted to avoid duplicates)
    unique_pairs = set(tuple(sorted(pair)) for pair in zip(df_filtered[col1_name], df_filtered[col2_name]))

    # Save unique pairs to a file
    with open(output_file, "w") as f:
        for pair in unique_pairs:
            f.write(f"{pair[0]},{pair[1]}\n")

    return unique_pairs

def print_and_count_unique_pairs(unique_pairs):
    print(f"Total unique pairs: {len(unique_pairs)}")
    for pair in unique_pairs:
        print(pair)

def protein_lens(msa_input):
    for protein in msa_input['sequences']:
        print(len(protein['protein']['sequence']))

def convert_fasta_to_msa_json(input_file):
    """
    Converts a FASTA file into a JSON formatted for MSA.

    Args:
        input_file (str): Path to the input FASTA file.
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()

    msa_input = {
        "name": os.path.basename(input_file).split('.')[0].upper(),
        "modelSeeds": [1],
        "sequences": [],
        "dialect": "alphafold3",
        "version": 1
    }

    current_protein = None
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_protein:
                msa_input["sequences"].append(current_protein)

            # Extract only the Uniprot accession ID (e.g., Q8BG89)
            protein_id = line.split('|')[1] if '|' in line else line.strip('>')

            current_protein = {"protein": {"id": protein_id, "sequence": ""}}
        else:
            if current_protein:
                current_protein["protein"]["sequence"] += line

    if current_protein:
        msa_input["sequences"].append(current_protein)

    # with open(output_file, 'w') as f:
    #     json.dump(msa_input, f, indent=4, separators=(',', ': '))
    return msa_input

if __name__ == '__main__':

    csv_file = "data_from_Dina/PSD_residue_pairs.csv"
    col1 = "acc_a"
    col2 = "acc_b"

    # Step 1: Extract unique protein pairs from CSV
    output_file = "data_from_Dina/unique_pairs.csv"
    pairs = extract_unique_pairs(csv_file, col1, col2, output_file)

    # Step 2: Get unique protein IDs from all pairs
    unique_sequences = set(seq for pair in pairs for seq in pair)
    print(f"Total unique sequences: {len(unique_sequences)}")

    # Step 3: Convert FASTA to MSA JSON
    msa_dict = convert_fasta_to_msa_json("data_from_Dina/synapse_all_seqs.txt")

    # Step 4: Filter sequences not in unique pairs
    before = len(msa_dict["sequences"])
    msa_dict["sequences"] = [
        seq for seq in msa_dict["sequences"]
        if seq["protein"]["id"] in unique_sequences]

    after = len(msa_dict["sequences"])
    print(f"Filtered sequences: {before - after} removed, {after} kept")

    # Step 5: Save filtered MSA JSON
    with open("data_from_Dina/synapse_pre_msa.json", 'w') as f:
        json.dump(msa_dict, f, indent=4, separators=(',', ': '))

