#!/usr/bin/env python3
"""
Delete pair directories from msa_output that contain one or more chains
originating from a '_low_' subunit.

Expects the following layout in the provided directory:
    subunits_info.json
    chain_id_mapping.json
    msa_output/

Usage:
    python filter_msa_dirs.py <directory>

Example:
    python filter_msa_dirs.py /path/to/data/
"""

import json
import os
import shutil
import sys


def get_low_chains(subunit_info: dict, chain_id_mapping: dict) -> set:
    low_chains = set()

    for new_chain_key, mapping in chain_id_mapping.items():
        original_chain = mapping["chain_id"]
        start_res      = mapping["start"]

        for subunit_name, subunit in subunit_info.items():
            if (original_chain in subunit["chain_names"]
                    and subunit["start_res"] == start_res):
                if "_low_" in subunit_name:
                    low_chains.add(new_chain_key.lower())
                break

    return low_chains


def get_dirs_to_delete(msa_output_dir: str, low_chains: set) -> list:
    to_delete = []

    for entry in os.scandir(msa_output_dir):
        if not entry.is_dir():
            continue

        parts = entry.name.lower().split("_")

        if len(parts) != 2:
            continue

        chain1, chain2 = parts
        if chain1 in low_chains or chain2 in low_chains:
            to_delete.append(entry.path)

    return to_delete


def main():
    if len(sys.argv) != 2:
        print("Usage: python filter_msa_dirs.py <directory>")
        sys.exit(1)

    base_dir = sys.argv[1]

    subunit_info_path  = os.path.join(base_dir, "subunits_info.json")
    chain_mapping_path = os.path.join(base_dir, "chain_id_mapping.json")
    msa_output_dir     = os.path.join(base_dir, "msa_output")

    # --- Validate inputs ---
    for path in (subunit_info_path, chain_mapping_path):
        if not os.path.isfile(path):
            print(f"ERROR: File not found: {path}")
            sys.exit(1)

    if not os.path.isdir(msa_output_dir):
        print(f"ERROR: msa_output directory not found: {msa_output_dir}")
        sys.exit(1)

    # --- Load inputs ---
    with open(subunit_info_path) as f:
        subunit_info = json.load(f)

    with open(chain_mapping_path) as f:
        chain_id_mapping = json.load(f)

    # --- Determine low chains ---
    low_chains = get_low_chains(subunit_info, chain_id_mapping)
    print(f"Low chains identified: {sorted(low_chains)}")

    # --- Find directories to delete ---
    to_delete = get_dirs_to_delete(msa_output_dir, low_chains)

    if not to_delete:
        print("No matching pair directories found. Nothing to delete.")
        return

    print(f"\nDirectories to be deleted ({len(to_delete)}):")
    for path in sorted(to_delete):
        print(f"  {path}")

    # --- Confirm and delete ---
    confirm = input("\nProceed with deletion? [y/N]: ").strip().lower()
    if confirm != "y":
        print("Aborted. Nothing was deleted.")
        return

    deleted, failed = 0, 0
    for path in to_delete:
        try:
            shutil.rmtree(path)
            print(f"  Deleted: {path}")
            deleted += 1
        except Exception as e:
            print(f"  ERROR deleting {path}: {e}")
            failed += 1

    print(f"\nDone. Deleted: {deleted}, Failed: {failed}")


if __name__ == "__main__":
    main()