import os
import json
import string
import itertools

def rename_subunit_chains(input_file, output_dir, mapping_file_path):
    """
    Renames chain names in a subunits_info JSON file with short alphabetical labels
    (e.g., A, B, ..., Z, AA, AB, ..., ZZ) and saves each subunit as a separate JSON file.

    Args:
        input_file (str): Path to the input subunits_info JSON file.
        output_dir (str): Directory where the individual JSON files will be saved.
        mapping_file_path (str): Path where the reverse mapping JSON will be saved.
    """
    def generate_labels():
        """Generator that yields unique uppercase alphabetical labels."""
        alphabet = string.ascii_uppercase
        for letter in alphabet:
            yield letter
        for pair in itertools.product(alphabet, repeat=2):
            yield ''.join(pair)

    # Load the input JSON file
    with open(input_file, 'r') as f:
        subunits_data = json.load(f)

    os.makedirs(output_dir, exist_ok=True)
    label_generator = generate_labels()
    chain_map = {}

    # Rename chain names and save each subunit as a separate file
    for subunit_name, subunit_info in subunits_data.items():
        new_chain_names = []
        for chain_name in subunit_info["chain_names"]:
            if chain_name not in chain_map:
                chain_map[chain_name] = next(label_generator)
            new_chain_names.append(chain_map[chain_name])
        subunit_info["chain_names"] = new_chain_names

        # Save the individual subunit JSON file
        output_file = os.path.join(output_dir, f"{subunit_name}.json")
        with open(output_file, 'w') as f:
            json.dump(subunit_info, f, indent=4)
        print(f"✅ Saved {output_file}")

    # Save the mapping file
    with open(mapping_file_path, 'w') as f:
        json.dump(chain_map, f, indent=4)

    print(f"✅ Chain mapping saved to {mapping_file_path}")

    if __name__ == '__main__':
        if len(sys.argv) < 2 or len(sys.argv) > 3:
            print("Usage: script <dir_path> [input_file]")
            sys.exit(1)

        dir_path = os.path.abspath(sys.argv[1])
        input_file = os.path.join(dir_path, 'subunits_info.json') if len(sys.argv) == 2 else os.path.abspath(
            sys.argv[2])

        print(f"Directory path: {dir_path}")
        print(f"Input file: {input_file}")

        rename_subunit_chains(input_file, dir_path, os.path.join(dir_path, 'chain_id_mapping.json'))