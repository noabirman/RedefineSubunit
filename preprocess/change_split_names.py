import os
import json
import sys

def update_json_names(complexes_folder):
    """
    Update the 'name' field and rename JSON files in the 'msa_output' subfolders
    to match the first 'id' in the 'sequences' list.

    Args:
        complexes_folder (str): Path to the 'complexes' folder.
    """
    for root, _, files in os.walk(complexes_folder):
        if "msa_output" in root:  # Process only msa_output subfolders
            for file in files:
                if file.endswith(".json"):
                    file_path = os.path.join(root, file)
                    # Check if the file exists
                    if not os.path.exists(file_path):
                        print(f"Warning: File {file_path} does not exist. Skipping.")
                        continue
                    # Read the JSON file
                    print(f"Processing file: {file_path}")
                    with open(file_path, "r") as f:
                        data = json.load(f)

                    # Get the first 'id' in the 'sequences' list
                    if "sequences" in data and data["sequences"]:
                        first_id = data["sequences"][0]["protein"]["id"]

                        # Update the 'name' field
                        data["name"] = first_id

                        # Write the updated JSON back to the file
                        with open(file_path, "w") as f:
                            json.dump(data, f, indent=4)

                        # Rename the file to match the first 'id'
                        new_file_path = os.path.join(root, f"{first_id}.json")
                        os.rename(file_path, new_file_path)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: script.py <COMPLEXES_FOLDER>")
        sys.exit(1)

    complexes_folder = sys.argv[1]
    update_json_names(complexes_folder)