#!/usr/bin/env python3
"""
Maps phosphorylation positions from global protein coordinates directly to
combfold subchain local coordinates, in a single step.

Inputs:
  - phospho_ccd.json: global phospho positions keyed by accession (e.g. P15498)
  - combfold subunits_info.json: subchain info with accession in name, start_res in global coords
  - combfold chain_id_mapping.json: maps renamed chain IDs (AY, AA...) to original chain names

Output:
  - phospho_final.json: phospho positions keyed by renamed chain ID (AY, AA...), local 1-indexed

Usage:
  python3 map_phospho_to_subchains.py \
    --phospho-input phospho_ccd.json \
    --subunits-info /path/to/combfold/subunits_info.json \
    --chain-mapping /path/to/combfold/chain_id_mapping.json \
    --output phospho_final.json
"""

import json
import argparse


VALID_PTM_RESIDUES = {
    "SEP": "S",   # phosphoserine
    "TPO": "T",   # phosphothreonine
    "PTR": "Y",   # phosphotyrosine
}


def extract_accession(subunit_name):
    """Extract accession ID from subunit name, e.g. 'P15498_high_1' -> 'P15498'."""
    return subunit_name.split("_")[0]


def map_phospho_to_subchains(phospho_ccd, subunits_info, chain_id_mapping):
    """
    Map global phospho positions to local subchain positions in one step.

    Args:
        phospho_ccd: dict keyed by accession, values are lists of dicts with 'ccd' and 'position'
        subunits_info: dict keyed by subunit name (e.g. 'P15498_high_1'), with 'sequence',
                       'chain_names', and 'start_res' (in global protein coordinates)
        chain_id_mapping: dict keyed by renamed chain ID (e.g. 'AY'), with 'chain_id' and
                          'start'/'end' (in global protein coordinates)

    Returns:
        dict keyed by renamed chain ID, values are lists of phospho dicts with local positions
    """
    # Build a lookup: original_chain_id -> list of subunit info dicts
    # e.g. 'F' -> [{'accession': 'P15498', 'sequence': ..., 'start_res': 1}, ...]
    chain_to_subunits = {}
    for subunit_name, info in subunits_info.items():
        accession = extract_accession(subunit_name)
        original_chain = info["chain_names"][0]
        start_res = info["start_res"]
        seq = info["sequence"]
        end_res = start_res + len(seq) - 1

        if original_chain not in chain_to_subunits:
            chain_to_subunits[original_chain] = []
        chain_to_subunits[original_chain].append({
            "accession": accession,
            "sequence": seq,
            "start_res": start_res,
            "end_res": end_res,
        })

    result = {}
    skipped_invalid = []

    for renamed_id, mapping in chain_id_mapping.items():
        original_chain = mapping["chain_id"]
        chain_start = mapping["start"]
        chain_end = mapping["end"]

        # Find the subunit(s) for this original chain
        subunits = chain_to_subunits.get(original_chain, [])
        if not subunits:
            continue

        # Get accession from the subunit (all subunits of same chain share accession)
        accession = subunits[0]["accession"]

        # Get global phospho positions for this accession
        phospho_list = phospho_ccd.get(accession, [])
        if not phospho_list:
            continue

        # Find the subunit that corresponds to this renamed chain's coordinate range
        matching_subunit = None
        for su in subunits:
            if su["start_res"] == chain_start and su["end_res"] == chain_end:
                matching_subunit = su
                break

        if matching_subunit is None:
            # fallback: find subunit that overlaps with chain_start..chain_end
            for su in subunits:
                if su["start_res"] <= chain_end and su["end_res"] >= chain_start:
                    matching_subunit = su
                    break

        if matching_subunit is None:
            continue

        sequence = matching_subunit["sequence"]
        local_phospho = []

        for phospho in phospho_list:
            global_pos = phospho["position"]
            ccd = phospho["ccd"]

            # Check if this position falls within this subchain's global range
            if not (chain_start <= global_pos <= chain_end):
                continue

            # Convert to local 1-indexed position
            local_pos = global_pos - chain_start + 1

            # Validate: check residue at this position matches the expected PTM residue
            expected_residue = VALID_PTM_RESIDUES.get(ccd)
            actual_residue = sequence[local_pos - 1] if 1 <= local_pos <= len(sequence) else None

            if expected_residue and actual_residue != expected_residue:
                skipped_invalid.append({
                    "chain": renamed_id,
                    "global_pos": global_pos,
                    "local_pos": local_pos,
                    "ccd": ccd,
                    "expected": expected_residue,
                    "actual": actual_residue,
                })
                continue

            local_phospho.append({
                "ccd": ccd,
                "position": local_pos,
                "original_position": global_pos,
            })

        if local_phospho:
            result[renamed_id] = local_phospho

    # Report any skipped invalid PTMs
    if skipped_invalid:
        print(f"\nWARNING: Skipped {len(skipped_invalid)} invalid PTM assignments "
              f"(residue mismatch):") #probably indicates an isoform.
        for s in skipped_invalid:
            print(f"  Chain {s['chain']}: global pos {s['global_pos']} (local {s['local_pos']}) "
                  f"- {s['ccd']} expects '{s['expected']}' but found '{s['actual']}'")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Map phospho positions from global accession coordinates to subchain local coordinates."
    )
    parser.add_argument("--phospho-input", required=True,
                        help="Path to phospho_ccd.json (global positions keyed by accession)")
    parser.add_argument("--subunits-info", required=True,
                        help="Path to combfold subunits_info.json")
    parser.add_argument("--chain-mapping", required=True,
                        help="Path to combfold chain_id_mapping.json")
    parser.add_argument("--output", required=True,
                        help="Output path for phospho_final.json")
    args = parser.parse_args()

    with open(args.phospho_input) as f:
        phospho_ccd = json.load(f)

    with open(args.subunits_info) as f:
        subunits_info = json.load(f)

    with open(args.chain_mapping) as f:
        chain_id_mapping = json.load(f)

    print(f"Loaded {len(phospho_ccd)} accessions with phospho sites")
    print(f"Loaded {len(subunits_info)} subunits")
    print(f"Loaded {len(chain_id_mapping)} renamed chains")

    result = map_phospho_to_subchains(phospho_ccd, subunits_info, chain_id_mapping)

    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\nSuccessfully mapped phospho sites to {len(result)} subchains")
    print(f"Output saved to: {args.output}")


if __name__ == "__main__":
    main()