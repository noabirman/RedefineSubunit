import json
import string
import itertools
import os

def generate_labels():
    """Generate labels A-Z, then AA, AB, AC... up to ZZ if needed."""
    alphabet = string.ascii_uppercase
    for letter in alphabet:
        yield letter
    for pair in itertools.product(alphabet, repeat=2):
        yield ''.join(pair)


def rename_protein_ids(input_file, output_file, mapping_file_path):
    with open(input_file, 'r') as f:
        msa_data = json.load(f)

    protein_map = {}
    label_generator = generate_labels()

    for seq in msa_data.get("sequences", []):
        original_id = seq["protein"]["id"]
        if original_id not in protein_map:
            protein_map[original_id] = next(label_generator)
        seq["protein"]["id"] = protein_map[original_id]

    # with open(mapping_file_path, 'w') as file:
    #     json.dump(protein_map, file, indent=4)
    reversed_map = {
        new_id: original_id.split("_")[0]
        for original_id, new_id in protein_map.items()
    }

    with open(mapping_file_path, 'w') as file:
        json.dump(reversed_map, file, indent=4)

    with open(output_file, 'w') as f:
        json.dump(msa_data, f, indent=4, separators=(',', ': '))
