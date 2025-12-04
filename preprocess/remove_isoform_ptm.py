import json
import os
import glob


def remove_modification_from_chain_j(json_path, position_to_remove=222):
    """
    Remove a specific modification from chain J in AlphaFold3 input JSON files.

    Args:
        json_path (str): Path to the input JSON file (e.g., L_J.json)
        position_to_remove (int): Position of modification to remove (default: 222)
    """
    # Load the JSON file
    with open(json_path, "r") as f:
        msa_file = json.load(f)

    # Track if we made any changes
    modifications_removed = 0

    # Process each sequence entry
    for seq_entry in msa_file["sequences"]:
        protein = seq_entry["protein"]
        pid = protein["id"]

        # Check if this is chain J (you can adjust this condition based on how chain J is identified)
        # Assuming the protein ID contains 'J' or you can check by other means
        if "J" in pid or pid.endswith("J"):
            # Filter out the modification at position 222
            if "modifications" in protein:
                original_count = len(protein["modifications"])
                protein["modifications"] = [
                    mod for mod in protein["modifications"]
                    if mod.get("ptmPosition") != position_to_remove
                ]
                new_count = len(protein["modifications"])
                modifications_removed += (original_count - new_count)

    # Save the modified JSON back to the same file (or to a new file)
    with open(json_path, "w") as f:
        json.dump(msa_file, f, indent=2)

    if modifications_removed > 0:
        print(f"Removed {modifications_removed} modification(s) at position {position_to_remove} from {json_path}")
    else:
        print(f"No modifications at position {position_to_remove} found in {json_path}")

    return modifications_removed


def process_all_j_files(input_dir, position_to_remove=222):
    """
    Process all *_J.json and J_*.json files in the input directory.

    Args:
        input_dir (str): Directory containing the JSON files
        position_to_remove (int): Position of modification to remove
    """
    # Find all files matching pattern *_J.json and J_*.json
    pattern1 = os.path.join(input_dir, "*_J.json")
    pattern2 = os.path.join(input_dir, "J_*.json")
    json_files = glob.glob(pattern1) + glob.glob(pattern2)

    if not json_files:
        print(f"No *_J.json files found in {input_dir}")
        return

    print(f"Found {len(json_files)} file(s) to process")
    total_removed = 0

    for json_file in json_files:
        removed = remove_modification_from_chain_j(json_file, position_to_remove)
        total_removed += removed

    print(f"\nTotal modifications removed: {total_removed}")


if __name__ == "__main__":
    # Set your input directory (where the *_J.json files are located)
    INPUT_DIR = "/cs/labs/dina/noabirman/tcellsUniprots/AF3_COMPLEX_RESULTS/af_triv_input"

    # Position to remove (default: 222)
    POSITION_TO_REMOVE = 222

    # Process all *_J.json files
    process_all_j_files(INPUT_DIR, POSITION_TO_REMOVE)

    # Or process a single file:
    # remove_modification_from_chain_j("/path/to/your/A_J.json", 222)