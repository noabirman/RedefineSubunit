import os
import json
import string
import itertools
import sys

def rename_subunits_and_create_msa_input(input_file, dir_path, min_length):
    """
    Renames subunits from a subunits_info JSON file, creates msa_input files for subunits
    with sequence length >= min_length, msa_output files for shorter ones, and saves a mapping file.

    Args:
        input_file (str): Path to the input subunits_info JSON file.
        dir_path (str): Directory where the msa_inputs folder and mapping file will be created.
        min_length (int): Minimum sequence length required to generate an msa_input file.
    """
    # Load the input JSON file
    with open(input_file, 'r') as f:
        subunits_data = json.load(f)

    # Create msa_inputs and msa_output directories
    msa_inputs_dir = os.path.join(dir_path, "msa_inputs")
    msa_outputs_dir = os.path.join(dir_path, "msa_output")
    os.makedirs(msa_inputs_dir, exist_ok=True)
    os.makedirs(msa_outputs_dir, exist_ok=True)

    # Initialize mapping dictionary
    mapping = {}

    # Generate labels A, B, ..., Z, AA, AB, ...
    def generate_labels():
        alphabet = string.ascii_uppercase
        for letter in alphabet:
            yield letter
        for pair in itertools.product(alphabet, repeat=2):
            yield ''.join(pair)

    label_generator = generate_labels()

    # Process each subunit
    for subunit_name, subunit_info in subunits_data.items():
        new_label = next(label_generator)
        sequence = subunit_info["sequence"]
        start_res = subunit_info["start_res"]
        sequence_length = len(sequence)
        end_res = start_res + sequence_length - 1

        # Update the mapping
        mapping[new_label] = {
            "chain_id": subunit_info["chain_names"][0],
            "start": start_res,
            "end": end_res
        }

        # Short sequence → msa_output
        if sequence_length < min_length:
            msa_output = {
                "name": new_label,
                "modelSeeds": [1],
                "sequences": [
                    {
                        "protein": {
                            "id": new_label,
                            "sequence": sequence,
                            "unpairedMsa": "",
                            "pairedMsa": "",
                            "templates": []
                        }
                    }
                ],
                "dialect": "alphafold3",
                "version": 1
            }
            label_lower = new_label.lower()
            output_subdir = os.path.join(msa_outputs_dir, label_lower)
            os.makedirs(output_subdir, exist_ok=True)
            output_file = os.path.join(output_subdir, f"{label_lower}_data.json")
            with open(output_file, 'w') as f:
                json.dump(msa_output, f, indent=4)
            print(f"Short sequence. Saved msa_output: {output_file}")

        # Long enough sequence → msa_input
        else:
            msa_input = {
                "name": new_label,
                "modelSeeds": [1],
                "sequences": [
                    {
                        "protein": {
                            "id": new_label,
                            "sequence": sequence
                        }
                    }
                ],
                "dialect": "alphafold3",
                "version": 1
            }
            output_file = os.path.join(msa_inputs_dir, f"{new_label}.json")
            with open(output_file, 'w') as f:
                json.dump(msa_input, f, indent=4)
            print(f"Saved msa_input: {output_file}")

    # Save the mapping file
    mapping_file_path = os.path.join(dir_path, "chain_id_mapping.json")
    with open(mapping_file_path, 'w') as f:
        json.dump(mapping, f, indent=4)
    print(f"Mapping file saved to {mapping_file_path}")

if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 4:
        print("Usage: script <combfold_dir_path> [min_length] [input_file]")
        sys.exit(1)

    combfold_dir_path = os.path.abspath(sys.argv[1])
    input_file = os.path.join(combfold_dir_path, 'subunits_info.json') if len(sys.argv) < 4 else os.path.abspath(sys.argv[3])
    min_length = int(sys.argv[2]) if len(sys.argv) == 3 else 10

    print(f"Directory path: {combfold_dir_path}")
    print(f"Input file: {input_file}")
    print(f"Minimum sequence length for msa_input: {min_length}")

    rename_subunits_and_create_msa_input(input_file, combfold_dir_path, min_length)
