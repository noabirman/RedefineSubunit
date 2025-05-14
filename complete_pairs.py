import os
import json

def normalize_pair(pair):
    return tuple(sorted(part.lower() for part in pair))

def find_missing_pairs(msa_pairs_dir, af_pairs_dir):
    # Collect pairs from MSA files
    msa_pairs = {normalize_pair(file_name.replace('.json', '').split('_')) for file_name in os.listdir(msa_pairs_dir) if file_name.endswith('.json')}

    # Collect pairs from AF directories
    af_pairs = {normalize_pair(dir_name.split('_')) for dir_name in os.listdir(af_pairs_dir) if os.path.isdir(os.path.join(af_pairs_dir, dir_name))}

    # Find the missing pairs
    missing_pairs = list(msa_pairs - af_pairs)
    return missing_pairs

def find_shared_chains(subunits_info_path):
    with open(subunits_info_path, 'r') as f:
        subunits_info = json.load(f)

    shared_chains = {}
    subunit_names = list(subunits_info.keys())
    for i, name1 in enumerate(subunit_names):
        for name2 in subunit_names[i + 1:]:
            common_chains = set(chain.lower() for chain in subunits_info[name1]['chain_names']).intersection(
                chain.lower() for chain in subunits_info[name2]['chain_names']
            )
            if common_chains:
                # Use a string key instead of a tuple
                key = f"{name1},{name2}"
                shared_chains[key] = list(common_chains)

    return shared_chains

def main(dir_path):
    output = {}
    for sub_dir in os.listdir(dir_path):
        sub_dir_path = os.path.join(dir_path, sub_dir)
        if not os.path.isdir(sub_dir_path):
            continue

        print(f"Processing subdirectory: {sub_dir}")
        msa_pairs_dir = os.path.join(sub_dir_path, 'msa_pairs')
        af_pairs_dir = os.path.join(sub_dir_path, 'af_pairs')
        subunits_info_path = os.path.join(sub_dir_path, 'subunits_info.json')

        if not (os.path.exists(msa_pairs_dir) and os.path.exists(af_pairs_dir) and os.path.exists(subunits_info_path)):
            print(f"Skipping {sub_dir}: Required files or directories are missing.")
            continue

        # Find missing pairs
        missing_pairs = find_missing_pairs(msa_pairs_dir, af_pairs_dir)

        # Find shared chains
        shared_chains = find_shared_chains(subunits_info_path)

        # Store results
        output[sub_dir] = {
            "missing_pairs": missing_pairs,
            "shared_chains": shared_chains
        }

        # Print results
        print(f"Subdirectory: {sub_dir}")
        print(f"  Missing pairs: {missing_pairs}")
        print(f"  Shared chains: {shared_chains}")

    # Write output to JSON file
    output_file = os.path.join(dir_path, 'summary.json')
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=4)

    print(f"\nSummary written to {output_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py <dir_path>")
    else:
        main(sys.argv[1])