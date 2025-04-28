import os
import json

def update_json_names(complexes_folder):
    """
    Update the 'name' field and rename JSON files in the 'msa_output' subfolders
    to match the first 'id' in the 'sequences' list.

    Args:
        complexes_folder (str): Path to the 'complexes' folder.
    """
    for root, _, files in os.walk(complexes_folder):
        print(root)
        if "msa_output" in root:  # Process only msa_output subfolders
            print(f"Processing folder: {root}")
            for file in files:
                print(file)
                if file.endswith(".json"):
                    file_path = os.path.join(root, file)

                    # Read the JSON file
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
                        print(f"Updated and renamed: {file_path} -> {new_file_path}")

if __name__ == "__main__":
    complexes_folder = "complexes"  # Replace with the path to your 'complexes' folder
    update_json_names(complexes_folder)