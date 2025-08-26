import json


#!/usr/bin/env python3
"""
af3_json_input_to_fasta.py

Convert a JSON file containing protein sequences into a FASTA file.

Usage:
    python af3_json_input_to_fasta.py input.json output.fasta

The JSON is expected to have the following structure:
{
    "sequences": [
        {
            "protein": {
                "id": "ProteinID",
                "sequence": "SEQUENCE"
            }
        },
        ...
    ]
}
"""

import sys
import os
import json


def json_to_fasta(input_file: str, output_file: str) -> None:
    """
    Convert a JSON file with protein sequences into FASTA format.

    Args:
        input_file (str): Path to the input JSON file.
        output_file (str): Path where the FASTA file will be saved.

    Raises:
        FileNotFoundError: If the input JSON file does not exist.
        KeyError: If the JSON structure does not contain the expected fields.
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    with open(input_file) as f:
        data = json.load(f)

    if "sequences" not in data:
        raise KeyError("Invalid JSON: 'sequences' field is missing")

    with open(output_file, "w") as fasta:
        for entry in data["sequences"]:
            try:
                pid = entry["protein"]["id"]
                seq = entry["protein"]["sequence"]
            except KeyError as e:
                raise KeyError(f"Invalid JSON entry: missing {e}")
            fasta.write(f">{pid}\n{seq}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python af3_json_input_to_fasta.py <input.json> <output.fasta>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    print(f"Reading input JSON: {input_path}")
    print(f"Writing FASTA to: {output_path}")

    try:
        json_to_fasta(input_path, output_path)
        print("Conversion completed successfully.")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
