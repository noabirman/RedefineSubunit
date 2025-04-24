import json
import sys
import os

# Documentation:
# This script converts an MSA JSON file into a subunit JSON file with the specified format.
# Usage:
#   python msa_to_subunits.py <msa_input_file> <output_folder>
# Arguments:
#   <msa_input_file>: Path to the input MSA JSON file.
#   <output_folder>: Path to the folder where the output JSON file will be saved.
# Output:
#   A JSON file named 'subunits_info.json' will be saved in the specified output folder.

def convert_msa_to_subunits(msa_input_file, output_folder):
    """
    Converts an MSA JSON file into a subunit JSON file with the specified format.

    Args:
        msa_input_file (str): Path to the input MSA JSON file.
        output_folder (str): Path to the folder where the output JSON file will be saved.

    The input MSA JSON file should have the following structure:
    {
        "sequences": [
            {
                "protein": {
                    "id": "CHAIN_ID",
                    "sequence": "SEQUENCE"
                }
            }
        ]
    }

    The output subunit JSON file will have the following structure:
    {
        "A0": {"name": "A0", "chain_names": ["A"], "start_res": 1, "sequence": "SEQUENCE"},
        ...
    }
    """
    # Load the input MSA JSON file
    with open(msa_input_file, 'r') as f:
        msa_data = json.load(f)

    # Initialize the output data structure
    subunits = {}

    # Extract and format the subunits
    for seq in (msa_data.get("sequences", [])):
        chain_id = seq["protein"]["id"]
        sequence = seq["protein"]["sequence"]
        subunits[chain_id] = {
            "name": chain_id,
            "chain_names": [chain_id],
            "start_res": 1,
            "sequence": sequence
        }

    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Save the output JSON file
    output_file = os.path.join(output_folder, "subunits_info.json")
    with open(output_file, 'w') as f:
        json.dump(subunits, f, indent=4)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python msa_to_subunits.py <msa_input_file> <output_folder>")
        sys.exit(1)
    msa_input_file, output_folder = sys.argv[1], sys.argv[2]
    convert_msa_to_subunits(msa_input_file, output_folder)