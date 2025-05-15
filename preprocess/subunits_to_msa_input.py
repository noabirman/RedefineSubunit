import json
import sys
import os

# Documentation:
# This script converts a JSON file containing subunit information into a JSON file formatted for MSA (Multiple Sequence Alignment).
# Usage:
#   python subunits_to_msa_input.py <input_file> <output_dir>
# Arguments:
#   <input_file>: Path to the input JSON file containing subunit information.
#   <output_dir>: Path to the directory where the output JSON file will be saved.
# Output:
#   A JSON file containing the subunit sequences formatted for MSA will be saved in the specified output directory.

def convert_subunits_to_msa_input(input_file, output_dir):
    """
    Converts a JSON file containing subunit information into a JSON file formatted for MSA.

    Args:
        input_file (str): Path to the input JSON file.
        output_dir (str): Path to the directory where the output JSON file will be saved.

    The input JSON file should have the following structure:
    {
        "subunit1": {
            "chain_names": ["A"],
            "sequence": "SEQUENCE1"
        },
        "subunit2": {
            "chain_names": ["B"],
            "sequence": "SEQUENCE2"
        }
    }

    The output JSON file will contain:
    - A "name" field derived from the input file name.
    - A "sequences" list, where each sequence is represented as a dictionary with "id" and "sequence".
    - Additional metadata fields such as "modelSeeds", "dialect", and "version".
    """
    # Load the input JSON file
    with open(input_file, 'r') as f:
        data = json.load(f)

    # Initialize the output data structure
    msa_input = {
        "name": input_file.split('/')[-1].split('.')[0],
        "modelSeeds": [1],
        "sequences": [],
        "dialect": "alphafold3",
        "version": 1
    }

    # Extract and format the sequences
    for subunit in data.values():
        sequence_info = {
            "id": subunit["chain_names"][0].upper(),
            "sequence": subunit["sequence"]
        }
        msa_input["sequences"].append({"protein": sequence_info})

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    # Save the output JSON file
    # Get the filename without the directory
    filename = os.path.basename(input_file)

    # Join the filename with the output directory
    output_file = os.path.join(output_dir, filename)
    with open(output_file, 'w') as f:
        json.dump(msa_input, f, indent=4)

if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: script <dir_path> [input_file]")
        sys.exit(1)

    dir_path = os.path.abspath(sys.argv[1])
    input_file = os.path.join(dir_path, 'subunits_info.json') if len(sys.argv) == 2 else os.path.abspath(sys.argv[2])

    convert_subunits_to_msa_input(input_file, dir_path)
