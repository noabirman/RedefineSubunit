import os
import json

def normalize_pair(pair):
    """
    Normalize a pair of strings by:
      - converting both parts to lowercase
      - sorting them alphabetically
    This ensures pairs like ("A", "B") and ("B", "A") are treated the same.

    Args:
        pair (list or tuple): Two elements representing a pair.

    Returns:
        tuple: Normalized pair as a lowercase, alphabetically sorted tuple.
    """
    return tuple(sorted(part.lower() for part in pair))

def find_missing_pairs(msa_pairs_dir, af_pairs_dir):
    """
    Compare MSA pairs and AF pairs to find missing pairs.

    MSA pairs come from .json files in msa_pairs_dir.
    AF pairs come from subdirectory names in af_pairs_dir.

    A pair is considered "missing" if it exists in MSA pairs but not in AF pairs.

    Args:
        msa_pairs_dir (str): Path to the MSA pairs directory.
        af_pairs_dir (str): Path to the AF pairs directory.

    Returns:
        list: List of missing pairs (each as a tuple).
    """
    # Collect pairs from MSA files
    msa_pairs = {normalize_pair(file_name.replace('.json', '').split('_')) for file_name in os.listdir(msa_pairs_dir) if file_name.endswith('.json')}

    # Collect pairs from AF directories
    af_pairs = {normalize_pair(dir_name.split('_')) for dir_name in os.listdir(af_pairs_dir) if os.path.isdir(os.path.join(af_pairs_dir, dir_name))}

    # Find the missing pairs
    missing_pairs = list(msa_pairs - af_pairs)
    return missing_pairs

def find_shared_chains(subunits_info_path):
    """
    Find subunit pairs that share at least one chain name.

    Reads subunits_info.json and checks each pair of subunits
    for common chain names.

    Args:
        subunits_info_path (str): Path to subunits_info.json file.

    Returns:
        dict: Keys are "subunit1,subunit2", values are lists of shared chains.
    """
    with open(subunits_info_path, 'r') as f:
        subunits_info = json.load(f)

    shared_chains = {}
    subunit_names = list(subunits_info.keys())
    for i, name1 in enumerate(subunit_names):
        for name2 in subunit_names[i + 1:]:
            common_chains = set(chain for chain in subunits_info[name1]['chain_names']).intersection(
                chain for chain in subunits_info[name2]['chain_names']
            )
            if common_chains:
                # Use a string key instead of a tuple
                key = f"{name1},{name2}"
                shared_chains[key] = list(common_chains)

    return shared_chains

def main(dir_path):
    """
    Process all subdirectories in dir_path to:
      1. Identify missing AF pairs compared to MSA pairs.
      2. Detect subunit pairs sharing chains.

    Results for each subdirectory are saved to summary.json in dir_path.

    Args:
        dir_path (str): Parent directory containing subdirectories to analyze.
    """
    output = {}
    subdirs_with_missing_pairs = []
    subdirs_with_shared_chains = []

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
        if missing_pairs:
            subdirs_with_missing_pairs.append(sub_dir)

        # Find shared chains
        shared_chains = find_shared_chains(subunits_info_path)
        if shared_chains:
            subdirs_with_shared_chains.append(sub_dir)

        # Store results
        output[sub_dir] = {
            "missing_pairs": missing_pairs,
            "shared_chains": shared_chains
        }

        # Print results
        print(f"Subdirectory: {sub_dir}")
        print(f"  Missing pairs: {missing_pairs}")
        print(f"  Shared chains: {shared_chains}")

    # Sort the output dictionary by subdirectory names
    sorted_output = dict(sorted(output.items()))

    # Add lists of subdirectories with missing pairs and shared chains
    sorted_output["subdirs_with_missing_pairs"] = subdirs_with_missing_pairs
    sorted_output["subdirs_with_shared_chains"] = subdirs_with_shared_chains

    # Write output to JSON file
    output_file = os.path.join(dir_path, 'summary.json')
    with open(output_file, 'w') as f:
        json.dump(sorted_output, f, indent=4)

    print(f"\nSummary written to {output_file}")