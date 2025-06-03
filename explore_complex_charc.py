import os
import json
from collections import defaultdict

def find_double_chain_proteins(base_dir="."):
    for entry in os.listdir(base_dir):
        path = os.path.join(base_dir, entry)
        if os.path.isdir(path):
            json_path = os.path.join(path, "subunits_info.json")
            if os.path.isfile(json_path):
                try:
                    with open(json_path, 'r') as f:
                        data = json.load(f)
                        for protein_info in data.values():
                            if isinstance(protein_info, dict):
                                chain_names = protein_info.get("chain_names", [])
                                if len(chain_names) > 1:
                                    print(f"{entry}: contain double chain")
                                    break  # Only print once per directory
                except Exception as e:
                    print(f"Failed to read {json_path}: {e}")
    print("finised find_double_chain_proteins")


def load_json(path):
    with open(path) as f:
        return json.load(f)


def analyze_complex(complex_path, complex_name):
    original_info = load_json(os.path.join(complex_path, "subunits_info.json"))
    redefine_info = load_json(os.path.join(complex_path, "combfold", "subunits_info.json"))
    mapping_file = os.path.join(complex_path, "chain_id_mapping.json")
    if not os.path.exists(mapping_file):
        return None  # Skip if no mapping

    chain_id_mapping = load_json(mapping_file)

    # Get original chain → sequence length
    original_chains = {}
    for entry in original_info.values():
        if not entry["chain_names"]:
            continue
        chain_id = entry["chain_names"][0]  # only take first
        original_chains[chain_id] = len(entry["sequence"])

    # Track which original chains were cut at preprocess
    preprocess_cuts = defaultdict(list)
    for new_chain, info in chain_id_mapping.items():
        orig_chain = info["chain_id"]
        preprocess_cuts[orig_chain].append((new_chain, info["start"], info["end"]))

    for val in preprocess_cuts.values():
        val.sort(key=lambda x: x[1])  # sort by start

    # Track redefine cuts
    redefine_cuts = defaultdict(list)
    for subunit in redefine_info.values():
        if not subunit["chain_names"]:
            continue
        orig_chain = subunit["chain_names"][0]  # again, first only
        redefine_cuts[orig_chain].append((subunit["name"], subunit["start_res"]))

    redefine_cut_summary = {
        c: sorted([(name, start) for name, start in v], key=lambda x: x[1])
        for c, v in redefine_cuts.items() if len(v) > 1
    }

    return {
        "complex_name": complex_name,
        "preprocess": dict(preprocess_cuts),
        "redefine": redefine_cut_summary,
        "num_redefine_chains": len(redefine_cut_summary)
    }


# def main(base_dir):
#     all_results = []
#     complex_names = []
#     redefine_summary = {}
#     total_redefined = 0
#
#     for complex_name in os.listdir(base_dir):
#         complex_path = os.path.join(base_dir, complex_name)
#         if not os.path.isdir(complex_path):
#             continue
#         result = analyze_complex(complex_path, complex_name)
#         if result:
#             all_results.append(result)
#             complex_names.append(complex_name)
#             if result["num_redefine_chains"] > 0:
#                 redefine_summary[complex_name] = result["num_redefine_chains"]
#                 total_redefined += 1
#
#     total_original_chains = 0
#     total_redefined_chains = 0
#
#     for r in all_results:
#         total_original_chains += len(r['preprocess'])  # count original chain_ids
#         total_redefined_chains += sum(len(v) for v in r['redefine'].values()) + \
#                                   (len(r['preprocess']) - len(r['redefine']))  # uncut chains are still present
#
#     print("===== GLOBAL SUMMARY =====")
#     print(f"Total complexes with redefine cuts: {total_redefined}")
#     print(f"Total original chains: {total_original_chains}")
#     print(f"Total chains after redefine: {total_redefined_chains}")
#     for cname, r in zip(complex_names, all_results):
#         total_chains = len(r['preprocess'])  # original chain count
#         redefined_chains = len(r['redefine'])  # how many chains were cut
#         if redefined_chains > 0:
#             print(f"  {cname}: {redefined_chains}/{total_chains} chain(s) redefined")
def main(base_dir):
    all_results = []
    complex_names = []
    redefine_summary = {}
    total_redefined = 0

    for complex_name in os.listdir(base_dir):
        complex_path = os.path.join(base_dir, complex_name)
        if not os.path.isdir(complex_path):
            continue
        result = analyze_complex(complex_path, complex_name)
        if result:
            all_results.append(result)
            complex_names.append(complex_name)
            if result["num_redefine_chains"] > 0:
                redefine_summary[complex_name] = result["num_redefine_chains"]
                total_redefined += 1

            # === PER-COMPLEX SUMMARY ===
            print(f"\n===== {complex_name} =====")
            print("Preprocess cuts:")
            for orig_chain, cuts in result['preprocess'].items():
                if len(cuts) > 1:
                    fragments = [f"{frag_id} ({start}-{end})" for frag_id, start, end in cuts]
                    print(f"  {orig_chain}: split into {len(cuts)} fragments -> {', '.join(fragments)}")

            print("Redefine cuts:")
            if result['redefine']:
                for orig_chain, cuts in result['redefine'].items():
                    fragments = [f"{frag_id} (start @ {start})" for frag_id, start in cuts]
                    print(f"  {orig_chain}: split into {len(cuts)} -> {', '.join(fragments)}")
            else:
                print("  No redefine cuts")

            # Chain count summary
            original_chain_count = len(result['preprocess'])
            redefine_chain_count = sum(len(v) for v in result['redefine'].values()) + \
                                   (original_chain_count - len(result['redefine']))
            print(f"Chains: {original_chain_count} → {redefine_chain_count}")

    # === GLOBAL SUMMARY ===
    total_original_chains = sum(len(r['preprocess']) for r in all_results)
    total_redefined_chains = sum(
        sum(len(v) for v in r['redefine'].values()) + (len(r['preprocess']) - len(r['redefine']))
        for r in all_results
    )

    print("\n===== GLOBAL SUMMARY =====")
    print(f"Total complexes with redefine cuts: {total_redefined}")
    print(f"Total original chains: {total_original_chains}")
    print(f"Total chains after redefine: {total_redefined_chains}")
    for cname, r in zip(complex_names, all_results):
        total_chains = len(r['preprocess'])  # original chain count
        redefined_chains = len(r['redefine'])  # how many were split
        if redefined_chains > 0:
            print(f"  {cname}: {redefined_chains}/{total_chains} chain(s) redefined")



if __name__ == "__main__":
    # Run the check from the current directory
    #find_double_chain_proteins(base_dir="/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2")
    main("/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2")


