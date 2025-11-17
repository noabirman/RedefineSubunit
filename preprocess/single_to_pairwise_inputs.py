#!/usr/bin/env python3
"""
Convert AlphaFold 3 single-chain MSA input JSONs into
pairwise MSA input JSONs, one per chain pair.
"""
import os
import json
import itertools
import argparse

def load_single_chain(path):
    """Extract chain ID, sequence, and modelSeeds from a single AF3 JSON."""
    with open(path) as f:
        data = json.load(f)

    seq = data["sequences"][0]["protein"]["sequence"]
    chain_id = data["sequences"][0]["protein"]["id"]
    model_seeds = data.get("modelSeeds", [1])

    return chain_id, seq, model_seeds


def main(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    json_files = [f for f in os.listdir(input_dir) if f.endswith(".json")]

    if not json_files:
        print("No JSON files found in input directory.")
        return

    print(f"Found {len(json_files)} input JSONs")

    # Load all chains
    chains = {}
    for f in json_files:
        path = os.path.join(input_dir, f)
        cid, seq, seeds = load_single_chain(path)
        chains[cid] = {"sequence": seq, "modelSeeds": seeds}

    chain_ids = sorted(chains.keys())
    print(f"Chains detected: {chain_ids}")

    # Build every pair
    for A, B in itertools.combinations(chain_ids, 2):

        out_path = os.path.join(output_dir, f"{A}_{B}.json")

        merged = {
            "name": f"{A}_{B}",
            "modelSeeds": chains[A]["modelSeeds"],   # preserve seed
            "sequences": [
                {"protein": {"id": A, "sequence": chains[A]["sequence"]}},
                {"protein": {"id": B, "sequence": chains[B]["sequence"]}}
            ],
            "dialect": "alphafold3",
            "version": 1
        }

        with open(out_path, "w") as f:
            json.dump(merged, f, indent=2)

        print(f"   Created {out_path}")

    print("\n Done! All pairwise MSA input JSONs written to:")
    print(f"  {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="folder with single-chain AF3 JSONs")
    parser.add_argument("output_dir", help="folder to write paired JSONs")
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
