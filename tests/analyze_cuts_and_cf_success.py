import os
import json
import csv
from collections import defaultdict
import pandas as pd
import re
from Bio.PDB import PDBParser



def load_json(path):
    with open(path) as f:
        return json.load(f)

def analyze_complex(complex_path, complex_name):
    subunit_info_path1 = os.path.join(complex_path, "subunits_info.json")
    subunit_info_path2 = os.path.join(complex_path, "combfold", "subunits_info.json")
    if not os.path.exists(subunit_info_path1) or not os.path.exists(subunit_info_path2):
        return None

    original_info = load_json(subunit_info_path1)
    redefine_info = load_json(subunit_info_path2)
    mapping_file = os.path.join(complex_path, "chain_id_mapping.json")
    if not os.path.exists(mapping_file):
        return None  # Skip if no mapping

    chain_id_mapping = load_json(mapping_file)

    # Get original chain -> sequence length
    original_chains = {}
    for entry in original_info.values():
        if not entry["chain_names"]:
            continue
        chain_id = entry["chain_names"][0]
        original_chains[chain_id] = len(entry["sequence"])

    # Track preprocess cuts
    preprocess_cuts = defaultdict(list)
    for new_chain, info in chain_id_mapping.items():
        orig_chain = info["chain_id"]
        preprocess_cuts[orig_chain].append((new_chain, info["start"], info["end"]))

    for val in preprocess_cuts.values():
        val.sort(key=lambda x: x[1])

    # Track redefine cuts
    redefine_cuts = defaultdict(list)
    for subunit in redefine_info.values():
        if not subunit["chain_names"]:
            continue
        orig_chain = subunit["chain_names"][0]
        redefine_cuts[orig_chain].append(
            (subunit["name"], subunit["start_res"], len(subunit["sequence"]))
        )

    # Compute chain counts
    original_chain_count = len(original_chains)
    redefine_chain_count = sum(len(v) for v in redefine_cuts.values())

    # Collect new subunit sizes (per chain)
    chain_fragments = []
    for chain, cuts in redefine_cuts.items():
        sizes = [size for (_, _, size) in sorted(cuts, key=lambda x: x[1])]
        if sizes:
            chain_fragments.append(f"{chain} -> {', '.join(str(s) for s in sizes)}")

    new_sizes = " | ".join(chain_fragments) if chain_fragments else "-"

    return {
        "complex_name": complex_name,
        "preprocess_cut": "Yes" if any(len(cuts) > 1 for cuts in preprocess_cuts.values()) else "No",
        "redefine_cut": "Yes" if any(len(cuts) > 1 for cuts in redefine_cuts.values()) else "No",
        "new_sizes": new_sizes,
        "chain_transition": f"{original_chain_count} -> {redefine_chain_count}",
    }

def check_cf_variants(root_dir):
    """Check success/failure for all combfold variants per complex"""
    variants = [
        "combfold",
        "combfold_all",
        "combfold_trivial",
        "combfold_us_trivial",
        "combfold_high"
    ]
    records = {}
    for complex_id in sorted(os.listdir(root_dir)):
        complex_dir = os.path.join(root_dir, complex_id)
        if not os.path.isdir(complex_dir):
            continue

        row = {}
        for variant in variants:
            variant_dir = os.path.join(complex_dir, variant, "results", "assembled_results")
            if os.path.exists(os.path.join(variant_dir, "output_clustered_0.pdb")):
                row[variant] = "success"
            else:
                row[variant] = "failed"
        records[complex_id] = row
    return records

def count_chains_in_subunits_info(subunits_info_path):
    with open(subunits_info_path) as f:
        data = json.load(f)
    total_chains = 0
    for subunit in data.values():
        chain_names = subunit.get("chain_names", [])
        total_chains += len(chain_names)
    return total_chains

def extract_cb_number(results_dir):
    if not os.path.isdir(results_dir):
        return None
    for fname in os.listdir(results_dir):
        match = re.match(r"cb_(\d+)_output\.res", fname)
        if match:
            return int(match.group(1))
    return None

def compute_chain_coverage(base_dir):
    coverage = {}
    for complex_name in os.listdir(base_dir):
        combfold_dir = os.path.join(base_dir, complex_name, "combfold")
        subunits_info_path = os.path.join(combfold_dir, "subunits_info.json")
        results_dir = os.path.join(combfold_dir, "results", "_unified_representation", "assembly_output")

        if not os.path.exists(subunits_info_path):
            continue

        total_chains = count_chains_in_subunits_info(subunits_info_path)
        cb_number = extract_cb_number(results_dir)

        if cb_number is None or total_chains == 0:
            coverage[complex_name] = None
        else:
            coverage[complex_name] = round(cb_number / total_chains, 3)  # fraction, 0–1
    return coverage

def count_expected_residues(subunits_info_path):
    with open(subunits_info_path) as f:
        data = json.load(f)
    total_residues = 0
    for subunit in data.values():
        sequence = subunit.get("sequence", "")
        chain_names = subunit.get("chain_names", [])
        total_residues += len(sequence) * len(chain_names)
    return total_residues

def find_cb_pdb_file(results_dir):
    if not os.path.isdir(results_dir):
        return None
    for fname in os.listdir(results_dir):
        if re.match(r"cb_\d+_output_0\.pdb", fname):
            return os.path.join(results_dir, fname)
    return None

def count_residues_in_pdb(pdb_path):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("model", pdb_path)
    except Exception as e:
        print(f" Failed to parse {pdb_path}: {e}")
        return None
    residue_count = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    residue_count += 1
    return residue_count

def compute_residue_coverage(base_dir):
    coverage = {}
    for complex_name in os.listdir(base_dir):
        combfold_dir = os.path.join(base_dir, complex_name, "combfold")
        subunits_info_path = os.path.join(combfold_dir, "subunits_info.json")
        pdb_dir = os.path.join(combfold_dir, "results", "assembled_results")

        if not os.path.exists(subunits_info_path):
            continue

        expected = count_expected_residues(subunits_info_path)
        pdb_path = find_cb_pdb_file(pdb_dir)
        if not pdb_path or expected == 0:
            coverage[complex_name] = None
            continue

        actual = count_residues_in_pdb(pdb_path)
        if actual is None:
            coverage[complex_name] = None
            continue

        coverage[complex_name] = round(actual / expected, 3)  # fraction 0–1
    return coverage

def main(base_dir, output_csv="complex_summary.csv"):
    # --- Gather combfold run status ---
    cf_status = check_cf_variants(base_dir)  # returns dict[complex] = {variant: success/failed}

    # --- Gather chain-wise coverage ---
    chain_cov_status = compute_chain_coverage(base_dir)  # dict[complex] = fraction

    # --- Gather residue coverage ---
    residue_cov_status = compute_residue_coverage(base_dir)  # dict[complex] = fraction

    # --- Column rename for short names ---
    col_rename = {
        "combfold": "us",
        "combfold_all": "all",
        "combfold_trivial": "trivial",
        "combfold_us_trivial": "us_trivial",
        "combfold_high": "high",
    }

    rows = []
    for complex_name in sorted(os.listdir(base_dir)):
        complex_path = os.path.join(base_dir, complex_name)
        if not os.path.isdir(complex_path):
            continue
        result = analyze_complex(complex_path, complex_name)
        if not result:
            continue

        # Merge combfold variant status
        if complex_name in cf_status:
            for old, new in col_rename.items():
                if old in cf_status[complex_name]:
                    result[new] = cf_status[complex_name][old]

        # Merge chain coverage
        if complex_name in chain_cov_status:
            cov = chain_cov_status[complex_name]
            result["chain_coverage"] = round(cov*100, 1) if cov is not None else "-"

        # Merge residue coverage
        if complex_name in residue_cov_status:
            cov = residue_cov_status[complex_name]
            result["residue_coverage"] = round(cov*100, 1) if cov is not None else "-"

        rows.append(result)

    # Headers
    extra_cols = list(col_rename.values()) + ["chain_coverage", "residue_coverage"]
    headers = ["complex_name", "preprocess_cut", "redefine_cut", "new_sizes", "chain_transition"] + extra_cols

    # --- Print table ---
    print(" | ".join(h.ljust(15) for h in headers))
    print("-" * 160)
    for r in rows:
        print(" | ".join(str(r.get(h, "-")).ljust(15) for h in headers))

    # --- Save CSV ---
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n Table saved to {output_csv}")



if __name__ == "__main__":
    root_dir = "/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2"
    main(root_dir, output_csv=os.path.join(root_dir, "complex_summary_full.csv"))
