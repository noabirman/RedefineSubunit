import os
import shutil
import json
import sys

def normalize_pair(pair):
    return tuple(sorted(pair))

def main(input_directory):
    msa_pairs_dir = os.path.join(input_directory, 'msa_pairs')
    af_pairs_dir = os.path.join(input_directory, 'af_pairs')
    json_file = os.path.join(input_directory, 'missing_pairs.json')
    
    if not os.path.exists(msa_pairs_dir) or not os.path.exists(af_pairs_dir) or not os.path.exists(json_file):
        print("Error: One or more required paths do not exist.")
        return
    
    # Read the JSON file
    with open(json_file, 'r') as f:
        missing_pairs = [normalize_pair(pair) for pair in json.load(f)]
    
    # Collect pairs from directories
    msa_pairs = {normalize_pair(dir_name.split('_')) for dir_name in os.listdir(msa_pairs_dir)}
    af_pairs = {normalize_pair(dir_name.split('_')) for dir_name in os.listdir(af_pairs_dir)}
    
    # Find the missing pairs
    real_missing_pairs = msa_pairs - af_pairs
    
    # Validation step
    if set(missing_pairs) != real_missing_pairs:
        print("Error: The pairs in the JSON do not match the actual missing pairs.")
        
        in_json_not_real = set(missing_pairs) - real_missing_pairs
        in_real_not_json = real_missing_pairs - set(missing_pairs)

        if in_json_not_real:
            print("\nPairs in JSON but not actually missing:")
            for pair in sorted(in_json_not_real):
                print(pair)
        
        if in_real_not_json:
            print("\nPairs actually missing but not listed in JSON:")
            for pair in sorted(in_real_not_json):
                print(pair)

        return
    
    # Create the output directory
    output_dir = os.path.join(input_directory, 'missing_msa_pairs')
    os.makedirs(output_dir, exist_ok=True)

    # Copy the missing pairs
    for pair in missing_pairs:
        dir_name = f"{pair[0]}_{pair[1]}"
        src_path = os.path.join(msa_pairs_dir, dir_name) if dir_name in os.listdir(msa_pairs_dir) else os.path.join(msa_pairs_dir, f"{pair[1]}_{pair[0]}")
        dst_path = os.path.join(output_dir, dir_name)
        shutil.copytree(src_path, dst_path)
    
    print(f"\nCopied {len(missing_pairs)} pairs to {output_dir}")
    
    # Read the JSON file
    with open(json_file, 'r') as f:
        missing_pairs = [normalize_pair(pair) for pair in json.load(f)]
    
    # Collect pairs from directories
    msa_pairs = {normalize_pair(dir_name.split('_')) for dir_name in os.listdir(msa_pairs_dir)}
    af_pairs = {normalize_pair(dir_name.split('_')) for dir_name in os.listdir(af_pairs_dir)}
    
    # Find the missing pairs
    real_missing_pairs = msa_pairs - af_pairs
    
    # Validation step
    if set(missing_pairs) != real_missing_pairs:
        print("Error: The pairs in the JSON do not match the actual missing pairs.")
        return
    
    # Create the output directory
    output_dir = os.path.join(input_directory, 'missing_msa_pairs')
    os.makedirs(output_dir, exist_ok=True)

    # Copy the missing pairs
    for pair in missing_pairs:
        dir_name = f"{pair[0]}_{pair[1]}"
        src_path = os.path.join(msa_pairs_dir, dir_name) if dir_name in os.listdir(msa_pairs_dir) else os.path.join(msa_pairs_dir, f"{pair[1]}_{pair[0]}")
        dst_path = os.path.join(output_dir, dir_name)
        shutil.copytree(src_path, dst_path)
    
    print(f"Copied {len(missing_pairs)} pairs to {output_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python sync_missing_pairs.py <input_directory>")
    else:
        main(sys.argv[1])

