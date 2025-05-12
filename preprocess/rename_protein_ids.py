import json
import string
import itertools
import os

def generate_labels():
    """
       Generator that yields unique uppercase alphabetical labels.

       Yields:
           str: A string label starting from 'A' to 'Z', then 'AA' to 'ZZ', etc.
       """
    alphabet = string.ascii_uppercase
    for letter in alphabet:
        yield letter
    for pair in itertools.product(alphabet, repeat=2):
        yield ''.join(pair)


def rename_protein_ids(input_file, output_file, mapping_file_path,iupred_mapping_path=None):
    """
        Renames protein IDs in an MSA JSON file with short alphabetical labels
        (e.g., A, B, ..., Z, AA, AB, ..., ZZ) and saves a reverse mapping.

        Args:
            input_file (str): Path to input MSA JSON file (must contain a "sequences" list with nested "protein" objects).
            output_file (str): Path where the modified MSA JSON will be written.
            mapping_file_path (str): Path where the reverse mapping JSON will be saved.
                                     Keys are short labels (e.g., 'A'), values are original IDs (up to the first underscore).

        Notes:
            - Only unique protein IDs are renamed.
            - The reverse mapping strips the underscore and anything after from original IDs.
        """
    with open(input_file, 'r') as f:
        msa_data = json.load(f)

    protein_map = {}
    label_generator = generate_labels()

    # Load optional start/end info from IUPred splitting
    iupred_info = {}
    if iupred_mapping_path:
        with open(iupred_mapping_path, 'r') as f:
            raw_map = json.load(f)
            for original_id, fragments in raw_map.items():
                for frag in fragments:
                    iupred_info[frag["id"]] = {
                        "chain_id": original_id,
                        "start": frag["start"],
                        "end": frag["end"]
                    }

    final_mapping = {}
    for seq in msa_data.get("sequences", []):
        original_id = seq["protein"]["id"]
        if original_id not in protein_map:
            protein_map[original_id] = next(label_generator)
        #seq["protein"]["id"] = protein_map[original_id]
        new_id = protein_map[original_id]
        seq["protein"]["id"] = new_id

        # Add start/end if available
        if original_id in iupred_info:
            final_mapping[new_id] = iupred_info[original_id]
        else:
            final_mapping[new_id] = {"chain_id": original_id.split("_")[0], "start": 1,
                                     "end": len(seq["protein"]["sequence"])}


# with open(mapping_file_path, 'w') as file:
#     #     json.dump(protein_map, file, indent=4)
#     reversed_map = {
#         new_id: original_id.split("_")[0]
#         for original_id, new_id in protein_map.items()
#     }

    with open(mapping_file_path, 'w') as file:
        json.dump(final_mapping, file, indent=4)

    with open(output_file, 'w') as f:
        json.dump(msa_data, f, indent=4, separators=(',', ': '))
