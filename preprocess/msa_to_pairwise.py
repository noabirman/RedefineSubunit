import json
import itertools
import os
import sys
from copy import deepcopy


def create_pairwise_msa(msa_folder: str, mapping_file: str, subunits_info_file: str, output_dir: str) -> None:
    """
    Create pairwise MSA files with self-pairs only for multi-chain originals
    """
    # Load mapping and subunit info
    with open(mapping_file, 'r') as f:
        mapping = json.load(f)

    with open(subunits_info_file, 'r') as f:
        subunits_info = json.load(f)

    # Get list of MSA files
    msa_files = [os.path.join(root, f) for root, _, files in os.walk(msa_folder) for f in files if f.endswith('.json')]
    subunits = {json.load(open(f))["sequences"][0]["protein"]["id"]: f for f in msa_files}    # Create an output directory
    os.makedirs(output_dir, exist_ok=True)

    # Generate combinations without self-pairs
    pairs = list(itertools.combinations(subunits, 2))


    # Add self-pairs for multi-chain originals
    for subunit in subunits:
        # Get original subunit name from mapping
        original_chain_name = mapping[subunit]
        # Find the original_name where subunits_info[key]['chain_names'] contains original_chain_name
        original_name = next(
            (key for key, info in subunits_info.items() if original_chain_name in info['chain_names']),
            None
        )
        if original_name and len(subunits_info.get(original_name, {}).get('chain_names', [])) > 1:
            print(f"Adding self-pair for subunit: {subunit}")
            pairs.append((subunit, subunit))

    # Create pair files
    for pair in pairs:
        sub1, sub2 = pair
        sub1_path = subunits[sub1]
        sub2_path = subunits[sub2]
        with open(sub1_path, 'r') as f1, open(sub2_path, 'r') as f2:
            seq1 = json.load(f1)
            seq2 = json.load(f2)

            # Initialize new dictionaries for seq_a and seq_b
            seq_a = {
                "protein": {
                    "id": "A",
                    "sequence": seq1['sequences'][0]['protein']['sequence'],
                    "unpairedMsa": seq1['sequences'][0]['protein'].get('unpairedMsa', ""),
                    "pairedMsa": "",
                    "templates": []
                }
            }
            seq_b = {
                "protein": {
                    "id": "B",
                    "sequence": seq2['sequences'][0]['protein']['sequence'],
                    "unpairedMsa": seq2['sequences'][0]['protein'].get('unpairedMsa', ""),
                    "pairedMsa": "",
                    "templates": []
                }
            }

        pair_name = f"{sub1}_{sub2}"
        output_path = os.path.join(output_dir, f"{pair_name}.json")

        with open(output_path, 'w') as f:
            json.dump({
                'name': pair_name,
                'modelSeeds': [1],
                'sequences': [seq_a, seq_b],
                'dialect': "alphafold3",
                'version': 1
            }, f, indent=4)

        print(f"Created {output_path}")


if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 4:
        print("Usage: script.py <MSA_FOLDER> [MAPPING_JSON] [SUBUNITS_INFO_JSON]")
        sys.exit(1)

    msa_folder = os.path.abspath(sys.argv[1])
    parent_dir = os.path.dirname(msa_folder)

    mapping_file = sys.argv[2] if len(sys.argv) > 2 else os.path.join(parent_dir, "chain_id_mapping.json")
    subunits_info_file = sys.argv[3] if len(sys.argv) > 3 else os.path.join(parent_dir, "subunits_info.json")
    output_dir = os.path.join(os.path.dirname(msa_folder), 'msa_pairs')

    create_pairwise_msa(msa_folder, mapping_file, subunits_info_file, output_dir)
