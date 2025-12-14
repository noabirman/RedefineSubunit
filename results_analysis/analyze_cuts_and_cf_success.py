import os
import json
import csv
from collections import defaultdict
import pandas as pd
import re
from Bio.PDB import PDBParser, MMCIFParser
# import warnings
# from Bio import BiopythonWarning
#warnings.simplefilter("ignore", BiopythonWarning)
import sys
COMPLEXES = ["8a3t", "8adl", "8cte", "8f5o", "7wff", "7e8t", "8hil", "7t3b", "7oba", "7uic", "7pkn", "7xho", "7zkq",
             "8a5o", "8fee", "8bel", "7qve", "7arc", "7ozn", "8adn", "7t2r", "7p2y", "7qru", "7use", "8e9g"]


def load_json(path):
    with open(path) as f:
        return json.load(f)

def total_sequence_length(subunits_info):
    """
    Sum the total sequence length across all chains in a subunits_info.json file.
    If a subunit has multiple chains, its sequence length is multiplied
    by the number of chains.

    Parameters
    ----------
    subunits_info_path : str
        Path to the subunits_info.json file.

    Returns
    -------
    int
        Total sequence length across all chains.
    """
    total_length = 0
    for subunit in subunits_info.values():
        seq = subunit.get("sequence", "")
        chain_names = subunit.get("chain_names", [])
        total_length += len(seq) * len(chain_names)

    return total_length

def analyze_complex(complex_path, complex_name, save_json=False):
    subunit_info_path1 = os.path.join(complex_path, "subunits_info.json")
    subunit_info_path2 = os.path.join(complex_path, "combfold", "subunits_info.json")
    if not os.path.exists(subunit_info_path1) or not os.path.exists(subunit_info_path2):
        return None

    original_info = load_json(subunit_info_path1)
    redefine_info = load_json(subunit_info_path2)
    mapping_file = os.path.join(complex_path, "chain_id_mapping.json")
    if not os.path.exists(mapping_file):
        return None  # Skip if no mapping

    orig_chain_num = count_chains_in_subunits_info(original_info)
    size_in_res = total_sequence_length(original_info)
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

    result = {
        "complex_name": complex_name,
        "preprocess_cut": "Yes" if any(len(cuts) > 1 for cuts in preprocess_cuts.values()) else "No",
        "redefine_cut": "Yes" if any(len(cuts) > 1 for cuts in redefine_cuts.values()) else "No",
        "new_sizes": new_sizes,
        "size": size_in_res,
        "number_of_chains": orig_chain_num,
        "subunit_transition": f"{len(original_info)} -> {len(redefine_info)}",
    }
    #  Optional save
    if save_json:
        out_path = os.path.join(complex_path, "analyze_cut.json")
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2)

    return result

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

def count_chains_in_subunits_info(subunits_info):
    total_chains = 0
    for subunit in subunits_info.values():
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

def count_chains(structure_file, structure_id="struct"):
    """
    Count the number of chains in a PDB or mmCIF structure file.

    Parameters
    ----------
    structure_file : str
        Path to the PDB or mmCIF file.
    structure_id : str, optional
        Identifier for the structure (default: "struct").

    Returns
    -------
    int
        Number of unique chains in the structure.
    """
    if structure_file.endswith(".pdb"):
        parser = PDBParser(QUIET=True)
    elif structure_file.endswith(".cif") or structure_file.endswith(".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError("File must be .pdb, .cif, or .mmcif")

    structure = parser.get_structure(structure_id, structure_file)

    chain_ids = set()
    for model in structure:
        for chain in model:
            chain_ids.add(chain.id)

    return len(chain_ids)

def compute_chain_coverage(base_dir):
    coverage = {}
    #for complex_name in os.listdir(base_dir):
    for complex_name in COMPLEXES:
        complex_dir = os.path.join(base_dir, complex_name)
        combfold_dir = os.path.join(complex_dir, "combfold")
        subunits_info_path = os.path.join(combfold_dir, "subunits_info.json")
        #results_dir = os.path.join(combfold_dir, "results", "_unified_representation", "assembly_output")
        pdb_dir = os.path.join(combfold_dir, "results", "assembled_results")
        pdb_path = find_models_path(complex_dir)[0]
        if not os.path.exists(subunits_info_path):
            continue
        subunits_info = load_json(subunits_info_path)
        if not os.path.exists(pdb_path):
            continue
        total_chains = count_chains_in_subunits_info(subunits_info)
        #cb_number = extract_cb_number(results_dir)
        num_chains_in_structure = count_chains(pdb_path)
        if total_chains == 0:
            coverage[complex_name] = None
        else:
            coverage[complex_name] = round(num_chains_in_structure / total_chains, 3)  # fraction, 0–1
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
    #for complex_name in os.listdir(base_dir):
    for complex_name in COMPLEXES:
        complex_dir = os.path.join(base_dir, complex_name)
        combfold_dir = os.path.join(complex_dir, "combfold")
        subunits_info_path = os.path.join(combfold_dir, "subunits_info.json")
        pdb_dir = os.path.join(combfold_dir, "results", "assembled_results")

        if not os.path.exists(subunits_info_path):
            continue

        expected = count_expected_residues(subunits_info_path)
        pdb_path = find_models_path(complex_dir)[0]
        if not pdb_path or expected == 0:
            coverage[complex_name] = None
            continue

        actual = count_residues_in_pdb(pdb_path)
        if actual is None:
            coverage[complex_name] = None
            continue

        coverage[complex_name] = round(actual / expected, 3)  # fraction 0–1
    return coverage

import os
import csv
import re
from disscussion import get_pdb_long_jumps

def analyze_structure_violations(structure_file_path):
    """
    Analyze all chains in a structure file for distance violations.
    Returns dict with chain violation counts and summary stats.
    """
    if not os.path.exists(structure_file_path):
        return {"total_violations": 0, "chains_with_violations": 0, "violation_details": "-"}

    # Determine file type and use appropriate parser
    file_ext = os.path.splitext(structure_file_path)[1].lower()
    if file_ext in ['.cif', '.mmcif']:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    try:
        structure = parser.get_structure('structure', structure_file_path)
        model = structure[0]
        chain_ids = [chain.id for chain in model]

        chain_violations = {}
        total_violations = 0

        for chain_id in chain_ids:
            violations_df = get_pdb_long_jumps(structure_file_path, chain_id)
            violations = len(violations_df)
            chain_violations[chain_id] = violations
            total_violations += violations

            # violations = get_pdb_long_jumps(structure_file_path, chain_id)
            # chain_violations[chain_id] = violations
            # total_violations += violations

        chains_with_violations = sum(1 for count in chain_violations.values() if count > 0)

        # Create violation details string
        if total_violations > 0:
            violation_details = ", ".join(
                [f"{chain}:{count}" for chain, count in chain_violations.items() if count > 0])
        else:
            violation_details = "None"

        return {
            "total_violations": total_violations,
            "chains_with_violations": chains_with_violations,
            "violation_details": violation_details,
            "total_chains": len(chain_ids)
        }

    except Exception as e:
        return {"total_violations": 0, "chains_with_violations": 0, "violation_details": "Error", "total_chains": 0}

def find_models_path(complex_dir):
    # returns the path for our structure model and combfold trivial structure model
    model_file = None
    cf_trivial_model_file = None
    structure_dir = os.path.join(complex_dir, "combfold", "results", "assembled_results")
    triv_structure_dir = os.path.join(complex_dir, "combfold_trivial", "results", "assembled_results")

    for filename in os.listdir(structure_dir):
        if filename.endswith("0.pdb"):
            model_file = os.path.join(structure_dir, filename)
            break
    for filename in os.listdir(triv_structure_dir):
        if filename.endswith("0.pdb"):
            cf_trivial_model_file = os.path.join(triv_structure_dir, filename)
            break
    return model_file, cf_trivial_model_file


def main(base_dir, output_csv="complex_summary.csv"):
    # --- Gather combfold run status ---
    cf_status = check_cf_variants(base_dir)  # returns dict[complex] = {variant: success/failed}

    # --- Gather chain-wise coverage ---
    chain_cov_status = compute_chain_coverage(base_dir)  # dict[complex] = fraction

    # --- Gather residue coverage ---
    residue_cov_status = compute_residue_coverage(base_dir)  # dict[complex] = fraction

    col_rename = {
        "combfold": "us",
        "combfold_all": "all",
        "combfold_trivial": "trivial",
        "combfold_us_trivial": "us_trivial",
        "combfold_high": "high",
    }

    rows = []
    #for complex_name in sorted(os.listdir(base_dir)):
    for complex_name in sorted(COMPLEXES):
        complex_path = os.path.join(base_dir, complex_name)
        if not os.path.isdir(complex_path):
            continue
        our_model, triv_model = find_models_path(complex_path)
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
            result["chain_coverage"] = round(cov * 100, 1) if cov is not None else "-"

        # Merge residue coverage
        if complex_name in residue_cov_status:
            cov = residue_cov_status[complex_name]
            result["residue_coverage"] = round(cov * 100, 1) if cov is not None else "-"

        # --- NEW: Add structure violation analysis ---
        # Structure file ends with "0.pdb" - find it in the complex directory

        violation_analysis = analyze_structure_violations(our_model)
        triv_violation_analysis = analyze_structure_violations(triv_model)

        # Add violation data to result
        result["distance_violations"] = violation_analysis["total_violations"]
        result["chains_with_violations"] = violation_analysis["chains_with_violations"]
        result["violation_details"] = violation_analysis["violation_details"]
        result["triv_distance_violations"] = triv_violation_analysis["total_violations"]
        result["triv_chains_with_violations"] = triv_violation_analysis["chains_with_violations"]
        result["triv_violation_details"] = triv_violation_analysis["violation_details"]

        rows.append(result)

    # Headers - added new columns for violations
    extra_cols = list(col_rename.values()) + ["chain_coverage", "residue_coverage"]
    violation_cols = ["distance_violations", "chains_with_violations", "violation_details","triv_distance_violations", "triv_chains_with_violations", "triv_violation_details"]
    headers = ["complex_name", "preprocess_cut", "redefine_cut", "new_sizes", "size",
               "number_of_chains", "subunit_transition"] + extra_cols + violation_cols

    # --- Print table ---
    print(" | ".join(h.ljust(20) for h in headers))
    print("-" * (21 * len(headers)))
    for r in rows:
        print(" | ".join(str(r.get(h, "-")).ljust(20) for h in headers))

    # --- Save CSV ---
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nTable saved to {output_csv}")

    # Print summary of violations for both structures
    total_complexes = len(rows)
    combfold_complexes_with_violations = sum(
        1 for r in rows if r.get("distance_violations", 0) not in [0, "-"])
    trivial_complexes_with_violations = sum(1 for r in rows if r.get("triv_distance_violations", 0) not in [0, "-"])

    print(f"\nViolation Summary:")
    print(
        f"  CombFold: {combfold_complexes_with_violations}/{total_complexes} complexes have distance threshold violations")
    print(
        f"  Trivial:  {trivial_complexes_with_violations}/{total_complexes} complexes have distance threshold violations")


if __name__ == "__main__":
    # root_dir = "/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2"
    # main(root_dir, output_csv=os.path.join(root_dir, "complex_summary_violations.csv"))
    #analyze_structure_violations("/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2/8a5o/combfold/results/assembled_results/cb_20_output_0.pdb")
    #analayze complex before combfold if you want after you need to take the first 3 rows
    complex_path = os.path.abspath(sys.argv[1])
    complex_name = os.path.abspath(sys.argv[2])
    analyze_complex(complex_path, complex_name, save_json=True)