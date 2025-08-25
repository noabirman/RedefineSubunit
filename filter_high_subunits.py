#!/usr/bin/env python3
"""
Filter a subunits_info.json file so that only entries whose key
contains '_high_' are kept.

Usage:
    python filter_high_subunits.py /path/to/subunits_info.json
"""

import sys
import os
import json


def load_json(path):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        sys.exit(f"File not found: {path}")
    except json.JSONDecodeError as e:
        sys.exit(f"Failed to parse JSON: {e}")


def save_json(data, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)


def keep_high_only(full_dict):
    """Return a new dict containing only keys with '_high_' (case-insensitive)."""
    return {k: v for k, v in full_dict.items() if "_high_" in k.lower()}


def main():
    if len(sys.argv) != 2:
        print("Usage: python filter_high_subunits.py /path/to/subunits_info.json")
        sys.exit(1)

    # Resolve the input path (handles ~ and makes it absolute)
    input_path = os.path.abspath(os.path.expanduser(sys.argv[1]))

    # --- read, filter, write --------------------------------------------------
    data = load_json(input_path)
    high_only = keep_high_only(data)

    if not high_only:
        sys.exit("No keys containing '_high_' were found - nothing to write.")

    output_path = os.path.join(os.path.dirname(input_path), "high_subunits_info.json")
    save_json(high_only, output_path)

    print(f"Saved {len(high_only)} high subunit(s) → {output_path}")


if __name__ == "__main__":
    main()
