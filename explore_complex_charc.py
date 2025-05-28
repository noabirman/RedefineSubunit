import os
import json

import os
import json

def find_double_chain_proteins(base_dir="."):
    for entry in os.listdir(base_dir):
        path = os.path.join(base_dir, entry)
        if os.path.isdir(path):
            json_path = os.path.join(path, "subunits_info.json")
            if os.path.isfile(json_path):
                try:
                    with open(json_path, 'r') as f:
                        data = json.load(f)
                        for protein_info in data.values():
                            if isinstance(protein_info, dict):
                                chain_names = protein_info.get("chain_names", [])
                                if len(chain_names) > 1:
                                    print(f"{entry}: contain double chain")
                                    break  # Only print once per directory
                except Exception as e:
                    print(f"Failed to read {json_path}: {e}")
    print("finised find_double_chain_proteins")


if __name__ == "__main__":
    # Run the check from the current directory
    find_double_chain_proteins(base_dir="/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2")
