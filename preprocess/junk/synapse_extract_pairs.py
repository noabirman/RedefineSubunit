from synapse_preprocess import extract_unique_pairs, print_and_count_unique_pairs

import json, os, shutil
import pandas as pd
from pathlib import Path
if __name__ == '__main__':
    # ---------- inputs ----------
    csv_file = "data_from_Dina/PSD_residue_pairs.csv"
    col1 = "acc_a"
    col2 = "acc_b"
    mapping_path    = "/cs/labs/dina/tsori/af3_example/synapseComplex/chain_id_mapping.json"
    json_source_dir = Path("/cs/labs/dina/tsori/af3_example/synapseComplex/msa_pairs")          # where the L_O.json etc. live
    json_target_dir = Path("/cs/labs/dina/tsori/af3_example/synapseComplex/filtered_msa_pairs")     # where you want to collect them
    json_target_dir.mkdir(exist_ok=True)

    # ---------- step 1 : get unique UniProt pairs ----------
    # Step 1: Extract unique protein pairs from CSV
    output_file = "data_from_Dina/unique_pairs_2.csv"
    uniprot_pairs = extract_unique_pairs(csv_file, col1, col2, output_file)

    # ---------- step 2 : load chain-to-UniProt map ----------
    with open(mapping_path) as f:
        chain2uni = json.load(f)

    # invert it → UniProt → list[chain]
    uni2chains = {}
    for chain, uni in chain2uni.items():
        uni2chains.setdefault(uni, []).append(chain)

    # ---------- step 3 : expand to all chain-letter pairs ----------
    chain_pairs = set()
    missing_raw  = []                       # keeps track of pairs for which neither file exists

    for u1, u2 in uniprot_pairs:
        for c1 in uni2chains.get(u1, []):
            for c2 in uni2chains.get(u2, []):
                p1 = f"{c1}_{c2}"           # e.g. L_O
                p2 = f"{c2}_{c1}"           # e.g. O_L

                if (json_source_dir / f"{p1}.json").exists():
                    chain_pairs.add(p1)
                elif (json_source_dir / f"{p2}.json").exists():
                    chain_pairs.add(p2)
                else:
                    missing_raw.append(f"{p1} / {p2}")

    print(f"{len(chain_pairs)} chain-letter pairs found.")
    if missing_raw:
        print("⚠️  No JSON for:", ", ".join(missing_raw))

    # ---------- step 4 : collect the corresponding JSON files ----------
    missing = []
    for pair in chain_pairs:
        src = json_source_dir / f"{pair}.json"
        dst = json_target_dir / src.name
        if src.exists():
            shutil.copy2(src, dst)
        else:
            missing.append(pair)

    if missing:
        print("⚠️  JSON not found for:", ", ".join(missing))
    else:
        print("✓ All pair files copied to", json_target_dir)

# # Load the CSV file of pairs (assuming format: col1,col2)
# csv_file = "data_from_Dina/unique_pairs.csv"
# df = pd.read_csv(csv_file, header=None, names=["chain1", "chain2"])
# #Invert mapping to go from protein_id to chain keys
# from collections import defaultdict
#
# protein_to_chains = defaultdict(list)
# for chain, prot_id in chain_id_mapping.items():
#     protein_to_chains[prot_id].append(chain)
#
# # Target directory for filtered JSON pairs
# json_target_dir = Path("filtered_pairs")
# json_target_dir.mkdir(exist_ok=True)
#
# for _, row in df.iterrows():
#     prot1 = row["chain1"]
#     prot2 = row["chain2"]
#
#     # Get all chains for each protein
#     chains1 = protein_to_chains.get(prot1, [])
#     chains2 = protein_to_chains.get(prot2, [])
#
#     # For each chain pair, create filenames and copy JSON files
#     for c1 in chains1:
#         for c2 in chains2:
#             pair_name = f"{c1}_{c2}"
#             json_file_name = f"{pair_name}.json"
#             source_path = Path("path_to_json_files") / json_file_name
#             target_path = json_target_dir / json_file_name
#
#             if source_path.exists():
#                 shutil.copy2(source_path, target_path)
#                 print(f"Copied {source_path} to {target_path}")
#             else:
#                 print(f"File not found: {source_path}")
