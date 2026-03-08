import json
import argparse


def first_run(phospho_in_orig_subunits, subunits_info, chain_id_mapping):
    """
    First preprocessing run: map phospho positions from original subunits to renamed chains.

    Args:
        phospho_in_orig_subunits: Dict with accession IDs as keys and phospho sites as values
        subunits_info: Dict with accession IDs containing chain_names
        chain_id_mapping: Dict with renamed chain IDs and their mapping info

    Returns:
        Dict with renamed chain IDs as keys and adjusted phospho positions as values
    """
    # Map original chain → list of phospho dicts
    map_chain_to_phos = {}
    for acc, phos_list in phospho_in_orig_subunits.items():
        chain_name = subunits_info[acc]["chain_names"][0]
        map_chain_to_phos[chain_name] = phos_list

    # Compute phospho positions after preprocessing
    phos_pos_after_pp = {}

    for renamed_chain_id, vals in chain_id_mapping.items():
        org_chain_id = vals["chain_id"]
        start = vals["start"]
        end = vals["end"]

        # all PTMs for original chain (list of dicts)
        all_phosp_in_org_chain = map_chain_to_phos.get(org_chain_id, [])

        # keep only PTMs within this slice
        relevant_phos = [
            p for p in all_phosp_in_org_chain
            if start <= p["position"] <= end
        ]

        # shift positions so they start at 1
        shifted = [
            {
                "ccd": p["ccd"],
                "position": p["position"] - start + 1
            }
            for p in relevant_phos
        ]

        phos_pos_after_pp[renamed_chain_id] = shifted

    return phos_pos_after_pp


def second_run(phospho_positions, subchain_mapping):
    """
    Second preprocessing run: map phospho positions from chains to subchains.

    Args:
        phospho_positions: Dict with chain IDs as keys and list of phospho sites as values
        subchain_mapping: Dict with subchain names as keys and chain info as values

    Returns:
        Dict with subchain names as keys and adjusted phospho positions as values
    """
    new_phospho_positions = {}

    # Iterate through each subchain
    for subchain_name, subchain_info in subchain_mapping.items():
        original_chain = subchain_info['chain_id']
        start = subchain_info['start']
        end = subchain_info['end']

        # Check if original chain has phospho positions
        if original_chain in phospho_positions:
            phospho_list = []

            # Check each phospho position in the original chain
            for phospho in phospho_positions[original_chain]:
                position = phospho['position']

                # Check if position falls within subchain range
                if start <= position <= end:
                    # Adjust position relative to subchain (1-indexed)
                    adjusted_position = position - start + 1

                    phospho_list.append({
                        'ccd': phospho['ccd'],
                        'position': adjusted_position,
                        'original_position': position
                    })

            # Only add subchain if it has phosphorylation sites
            if phospho_list:
                new_phospho_positions[subchain_name] = phospho_list

    return new_phospho_positions


def main():
    # first run: data/phospho_ccd.json, data/subunits_info.json, data/chain_id_mapping.json
    # second run: phospho_after_first_preprocess.json, combfold/chain_id_mapping.json
    parser = argparse.ArgumentParser(
        description='Map phosphorylation positions through preprocessing pipeline'
    )
    parser.add_argument(
        '--mode',
        type=int,
        choices=[1, 2],
        required=True,
        help='Mode 1: first preprocessing run, Mode 2: second preprocessing run'
    )
    parser.add_argument(
        '--phospho-input',
        type=str,
        required=True,
        help='Input phospho JSON file (mode 1: data/phospho_ccd.json, mode 2: phospho_after_first_preprocess.json)'
    )
    parser.add_argument(
        '--mapping-input',
        type=str,
        required=True,
        help='Chain mapping JSON file (mode 1: data/chain_id_mapping.json, mode 2: combfold/chain_id_mapping.json)'
    )
    parser.add_argument(
        '--subunits-input',
        type=str,
        default=None,
        help='Subunits info JSON file (required for mode 1: data/subunits_info.json)'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output JSON file path'
    )

    args = parser.parse_args()

    # Load phospho input
    with open(args.phospho_input, "r") as f:
        phospho_data = json.load(f)

    # Load chain mapping
    with open(args.mapping_input, "r") as f:
        chain_mapping = json.load(f)

    # Process based on mode
    if args.mode == 1:
        if not args.subunits_input:
            parser.error("--subunits-input is required for mode 1")

        # Load subunits info
        with open(args.subunits_input, "r") as f:
            subunits_info = json.load(f)

        print("Running first preprocessing...")
        result = first_run(phospho_data, subunits_info, chain_mapping)
    else:  # mode == 2
        print("Running second preprocessing...")
        result = second_run(phospho_data, chain_mapping)

    # Save result
    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)

    print(f"Successfully processed {len(result)} chains/subchains with phosphorylation sites")
    print(f"Output saved to: {args.output}")


if __name__ == '__main__':
    main()
    #run the second with 2 data2/phospho_after_first_preprocess.json data2/chain_id_mapping.json