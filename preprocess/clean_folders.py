import os
import re
import sys
import shutil
from pathlib import Path
from datetime import datetime
import json


def remove_duplicates(input_dir, output_dir):
    # Match folders like a_b or a_b_YYYYMMDD_HHMMSS
    folder_pattern = re.compile(r"^([a-zA-Z]+_[a-zA-Z]+)(?:_(\d{8}_\d{6}))?$")

    groups = {}
    print("\n🔍 Scanning folders...")
    for entry in input_dir.iterdir():
        if entry.is_dir():
            match = folder_pattern.match(entry.name)
            if match:
                base = match.group(1)
                timestamp_str = match.group(2)
                timestamp = None
                if timestamp_str:
                    try:
                        timestamp = datetime.strptime(timestamp_str, "%Y%m%d_%H%M%S")
                    except ValueError:
                        print(f"  ⚠️ Skipping invalid timestamp: {entry.name}")
                        continue
                else:
                    # No timestamp: treat as oldest possible
                    timestamp = datetime.min

                print(f"  ✔ Matched: {entry.name} (base='{base}', timestamp='{timestamp}')")
                groups.setdefault(base, []).append((entry, timestamp))
            else:
                print(f"  ⚠️ Skipping non-matching folder: {entry.name}")

    print("\n📊 Processing groups to keep newest folder only...")
    for base, items in groups.items():
        print(f"\n➡️ Group: {base}")
        for entry, ts in items:
            print(f"   - {entry.name} (timestamp={ts})")

        # Sort by timestamp, newest last
        items.sort(key=lambda x: x[1])
        to_keep = items[-1][0]
        to_move = [item[0] for item in items if item[0] != to_keep]

        print(f"   ✅ Keeping: {to_keep.name}")
        for folder in to_move:
            print(f"   🗂️ Moving: {folder.name} → {output_dir}")
            shutil.move(str(folder), output_dir)
        new_name = input_dir / base
        if to_keep != new_name:
            print(f"   📝 Renaming: {to_keep.name} → {new_name.name}")
            to_keep.rename(new_name)
            to_keep = new_name

    print("\n✅ Done.")


def process_pairs(input_dir, output_dir):
    # Collect folder names
    folder_names = {entry.name for entry in input_dir.iterdir() if entry.is_dir()}

    # Check for pairs like a_b and b_a
    processed_pairs = set()
    for name in folder_names:
        if name in processed_pairs:
            continue
        parts = name.split('_')
        if len(parts) == 2:
            reverse_name = f"{parts[1]}_{parts[0]}"
            if reverse_name in folder_names:
                # Alphabetically determine which to move
                to_keep, to_move = sorted([reverse_name, name])
                print(f"🔄 Found pair: {name} and {reverse_name}")
                print(f"   ✅ Keeping: {to_keep}")
                print(f"   🗂️ Moving: {to_move} → {output_dir}")
                shutil.move(str(input_dir / to_move), output_dir)
                processed_pairs.add(name)
                processed_pairs.add(reverse_name)


def check_missing_pairs(input_dir, chain_names):
    # Get existing folder names
    folder_names = {entry.name for entry in input_dir.iterdir() if entry.is_dir()}

    # Generate all possible pairs
    missing_pairs = []
    for i, chain1 in enumerate(chain_names):
        for chain2 in chain_names[i + 1:]:
            pair1 = f"{chain1}_{chain2}"
            pair2 = f"{chain2}_{chain1}"

            if pair1 not in folder_names and pair2 not in folder_names:
                missing_pairs.append((chain1, chain2))

    if missing_pairs:
        print("\n❌ Missing pairs:")
        for chain1, chain2 in missing_pairs:
            print(f"   - {chain1}_{chain2} or {chain2}_{chain1}")
    else:
        print("\n✅ All possible pairs exist")


def load_chain_names(mapping_path):
    try:
        with open(mapping_path, 'r') as f:
            mapping = json.load(f)
            return [key.lower() for key in mapping.keys()]
    except FileNotFoundError:
        print(f"❌ Error: Mapping file not found at {mapping_path}")
        return []
    except json.JSONDecodeError:
        print(f"❌ Error: Invalid JSON in mapping file {mapping_path}")
        return []


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python clean_folders.py /path/to/input [/path/to/mapping.json]")
        sys.exit(1)

    input_path = Path(sys.argv[1]).resolve()
    if not input_path.is_dir():
        print(f"❌ Error: {input_path} is not a directory.")
        sys.exit(1)

    # Default mapping path or user-specified
    mapping_path = Path(sys.argv[2]) if len(sys.argv) > 2 else input_path.parent / 'chain_id_mapping.json'

    print(f"📂 Input directory: {input_path}")
    print(f"📄 Mapping file: {mapping_path}")
    print("Number of folders before: ", len(os.listdir(input_path)))

    # Output folder: one level above input, named "doubles"
    output_path = input_path.parent / "doubles"
    output_path.mkdir(exist_ok=True)
    print(f"📁 Archive directory: {output_path}")

    # Process folders
    remove_duplicates(input_path, output_path)
    process_pairs(input_path, output_path)

    # Check missing pairs using chain names from mapping
    chain_names = load_chain_names(mapping_path)
    if chain_names:
        print(f"\n🔍 Checking missing pairs for chains: {', '.join(chain_names)}")
        check_missing_pairs(input_path, chain_names)
    else:
        print("\n⚠️ No chain names loaded from mapping file")

    print("Number of folders after: ", len(os.listdir(input_path)))


