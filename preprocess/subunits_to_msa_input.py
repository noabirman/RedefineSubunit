import json
import sys
import os

def convert_subunits_to_msa_input(input_file, output_dir):
    """
    Converts a JSON file containing subunit information into a JSON file formatted for MSA.

    Args:
        input_file (str): Path to the input JSON file.
        output_dir (str): Path to the directory where the output JSON file will be saved.
    """
    print(f"Loading input file: {input_file}")
    with open(input_file, 'r') as f:
        data = json.load(f)

    print("Initializing MSA input structure...")
    msa_input = {
        "name": os.path.basename(input_file).split('.')[0],
        "modelSeeds": [1],
        "sequences": [],
        "dialect": "alphafold3",
        "version": 1
    }

    print("Extracting and formatting sequences...")
    for subunit_name, subunit in data.items():
        print(f"Processing subunit: {subunit_name}")
        sequence_info = {
            "id": subunit["chain_names"][0].upper(),
            "sequence": subunit["sequence"]
        }
        msa_input["sequences"].append({"protein": sequence_info})

    print(f"Creating output directory if it doesn't exist: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, os.path.basename(input_file))
    print(f"Saving MSA input to: {output_file}")
    with open(output_file, 'w') as f:
        json.dump(msa_input, f, indent=4)

    print("Conversion completed successfully.")

if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: script <dir_path> [input_file]")
        sys.exit(1)

    dir_path = os.path.abspath(sys.argv[1])
    input_file = os.path.join(dir_path, 'subunits_info.json') if len(sys.argv) == 2 else os.path.abspath(sys.argv[2])

    print(f"Directory path: {dir_path}")
    print(f"Input file: {input_file}")

    convert_subunits_to_msa_input(input_file, dir_path)