import os
import json
import string
import itertools

def rename_subunits_and_create_msa_input(input_file, dir_path):
    """
    Renames subunits in a subunits_info JSON file to A, B, C, ..., Z, AA, AB, etc.,
    creates a mapping file, and saves renamed subunits as individual msa_input JSON files
    in a new folder called msa_inputs.

    Args:
        input_file (str): Path to the input subunits_info JSON file.
        dir_path (str): Directory where the msa_inputs folder and mapping file will be created.
    """
    # Load the input JSON file
    with open(input_file, 'r') as f:
        subunits_data = json.load(f)

    # Create the msa_inputs directory
    msa_inputs_dir = os.path.join(dir_path, "msa_inputs")
    os.makedirs(msa_inputs_dir, exist_ok=True)

    # Initialize variables
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
        # Generate the new label
        new_label = next(label_generator)

        # Calculate the end residue
        start_res = subunit_info["start_res"]
        sequence_length = len(subunit_info["sequence"])
        end_res = start_res + sequence_length - 1

        # Update the mapping
        mapping[new_label] = {
            "chain_id": subunit_info["chain_names"][0],
            "start": start_res,
            "end": end_res
        }

        # Create the msa_input JSON structure
        msa_input = {
            "name": new_label,
            "modelSeeds": [1],
            "sequences": [
                {
                    "protein": {
                        "id": new_label,
                        "sequence": subunit_info["sequence"]
                    }
                }
            ],
            "dialect": "alphafold3",
            "version": 1
        }

        # Save the msa_input JSON file
        output_file = os.path.join(msa_inputs_dir, f"{new_label}.json")
        with open(output_file, 'w') as f:
            json.dump(msa_input, f, indent=4)
        print(f"✅ Saved {output_file}")

    # Save the mapping file
    mapping_file_path = os.path.join(dir_path, "chain_id_mapping.json")
    with open(mapping_file_path, 'w') as f:
        json.dump(mapping, f, indent=4)
    print(f"✅ Mapping file saved to {mapping_file_path}")

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: script <dir_path> [input_file]")
        sys.exit(1)

    dir_path = os.path.abspath(sys.argv[1])
    input_file = os.path.join(dir_path, 'subunits_info.json') if len(sys.argv) == 2 else os.path.abspath(sys.argv[2])

    print(f"Directory path: {dir_path}")
    print(f"Input file: {input_file}")

    rename_subunits_and_create_msa_input(input_file, dir_path)