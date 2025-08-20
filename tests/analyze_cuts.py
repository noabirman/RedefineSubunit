import os
import json
import csv
from collections import defaultdict

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

    # Track preprocess cuts
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
        orig_chain = subunit["chain_names"][0]
        redefine_cuts[orig_chain].append((subunit["name"], subunit["start_res"], len(subunit["sequence"])))

    redefine_cut_summary = {
        c: sorted(v, key=lambda x: x[1])
        for c, v in redefine_cuts.items() if len(v) > 1
    }

    # === Compute chain counts ===
    original_chain_count = len(original_chains)
    redefine_chain_count = sum(len(v) for v in redefine_cuts.values())

    # === Collect new subunit sizes (per chain) ===
    chain_fragments = []
    for chain, cuts in redefine_cuts.items():
        sizes = [size for (_, _, size) in sorted(cuts, key=lambda x: x[1])]
        if sizes:
            chain_fragments.append(f"{chain} → {', '.join(str(s) for s in sizes)}")

    new_sizes = " | ".join(chain_fragments) if chain_fragments else "-"

    return {
        "complex_name": complex_name,
        "preprocess_cut": "Yes" if any(len(cuts) > 1 for cuts in preprocess_cuts.values()) else "No",
        "redefine_cut": "Yes" if any(len(cuts) > 1 for cuts in redefine_cuts.values()) else "No",
        "new_sizes": new_sizes,
        "chain_transition": f"{original_chain_count} → {redefine_chain_count}",
    }

def main(base_dir, output_csv="complex_summary.csv"):
    rows = []
    for complex_name in sorted(os.listdir(base_dir)):
        complex_path = os.path.join(base_dir, complex_name)
        if not os.path.isdir(complex_path):
            continue
        result = analyze_complex(complex_path, complex_name)
        if result:
            rows.append(result)

    # === Print as table ===
    print(f"{'Complex':<15} {'Preprocess':<12} {'Redefine':<10} {'New Sizes':<60} {'Chains':<10}")
    print("-"*120)
    for r in rows:
        print(f"{r['complex_name']:<15} {r['preprocess_cut']:<12} {r['redefine_cut']:<10} "
              f"{r['new_sizes']:<60} {r['chain_transition']:<10}")

    # === Save to CSV ===
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["complex_name", "preprocess_cut", "redefine_cut", "new_sizes", "chain_transition"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n✅ Table saved to {output_csv}")

if __name__ == "__main__":
    main("/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2")
