import os
import json
from collections import defaultdict
import re
from Bio.PDB import PDBParser
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

def missing_chains_in_cf(base_dir):
    print("===== MISSING CHAIN SUMMARY =====")
    for complex_name in os.listdir(base_dir):
        combfold_dir = os.path.join(base_dir, complex_name, "combfold")
        subunits_info_path = os.path.join(combfold_dir, "subunits_info.json")
        results_dir = os.path.join(combfold_dir, "results", "_unified_representation", "assembly_output")

        if not os.path.exists(subunits_info_path):
            continue

        total_chains = count_chains_in_subunits_info(subunits_info_path)
        cb_number = extract_cb_number(results_dir)

        if cb_number is None:
            print(f"{complex_name}: ❌ No cb_*_output.res file found")
            continue

        missing = total_chains - cb_number
        print(f"{complex_name}: subunits_info = {total_chains}, cb_output = {cb_number}, missing = {missing}, existing: {cb_number/total_chains}")

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
        print(f"❌ Failed to parse {pdb_path}: {e}")
        return None
    residue_count = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    residue_count += 1
    return residue_count

def compare_residue_counts(base_dir):
    print("===== RESIDUE COUNT COMPARISON =====")
    for complex_name in os.listdir(base_dir):
        combfold_dir = os.path.join(base_dir, complex_name, "combfold")
        subunits_info_path = os.path.join(combfold_dir, "subunits_info.json")
        pdb_dir = os.path.join(combfold_dir, "results", "assembled_results")

        if not os.path.exists(subunits_info_path):
            continue

        expected_residues = count_expected_residues(subunits_info_path)
        pdb_path = find_cb_pdb_file(pdb_dir)
        if not pdb_path:
            print(f"{complex_name}: ❌ No cb_*_output_0.pdb found")
            continue

        actual_residues = count_residues_in_pdb(pdb_path)
        if actual_residues is None:
            continue

        missing = expected_residues - actual_residues
        percent_existing = 100 * actual_residues / expected_residues if expected_residues > 0 else 0

        print(f"{complex_name}: expected = {expected_residues}, actual = {actual_residues}, missing = {missing}, coverage: ({percent_existing:.1f}%)")

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
            if original_chain_count > 0:
                percent_increase = 100 * (redefine_chain_count - original_chain_count) / original_chain_count
                print(f"Chains: {original_chain_count} → {redefine_chain_count} (+{percent_increase:.1f}% increase)")
            else:
                print(f"Chains: {original_chain_count} → {redefine_chain_count} (N/A% increase)")

    # === GLOBAL SUMMARY ===
    total_original_chains = sum(len(r['preprocess']) for r in all_results)
    total_redefined_chains = sum(
        sum(len(v) for v in r['redefine'].values()) + (len(r['preprocess']) - len(r['redefine']))
        for r in all_results
    )

    print("\n===== GLOBAL SUMMARY =====")
    print(f"Total complexes with redefine cuts: {total_redefined}")
    print(f"Total original chains: {total_original_chains}")
    if total_original_chains > 0:
        percent_increase_global = 100 * (total_redefined_chains - total_original_chains) / total_original_chains
        print(
            f"Chain count increased from {total_original_chains} → {total_redefined_chains} (+{percent_increase_global:.1f}% increase)")
    else:
        print(f"Chain count increased from 0 → {total_redefined_chains} (N/A% increase)")

    for cname, r in zip(complex_names, all_results):
        total_chains = len(r['preprocess'])  # original chain count
        redefined_chains = len(r['redefine'])  # how many were split
        if redefined_chains > 0:
            print(f"  {cname}: {redefined_chains}/{total_chains} chain(s) redefined")
    missing_chains_in_cf(base_dir)
    compare_residue_counts(base_dir)

import matplotlib.pyplot as plt
# Count high and acceptable scores for "our" and "Ben's" TM-scores
def count_quality(scores, threshold):
    return sum(score > threshold for score in scores)

def plot_results_tm_score(json_path):
    with open(json_path) as f:
        data = json.load(f)
    our_scores = [entry["our_tm_score"] for entry in data]
    ben_scores = [entry["ben_best_tm_score"] for entry in data]

    total = len(data)

    our_scores = [entry["our_tm_score"] for entry in data]
    ben_scores = [entry["ben_best_tm_score"] for entry in data]
    total = len(data)

    # Calculate counts and percentages
    #our_high = count_quality(our_scores, 0.8)
    our_high = sum(d["our_tm_score"] > 0.8 for d in data)
    #ben_high = count_quality(ben_scores, 0.8)
    ben_high = sum(d["ben_best_tm_score"] > 0.8 for d in data)


    #our_acceptable = count_quality(our_scores, 0.7)
    #ben_acceptable = count_quality(ben_scores, 0.7)
    our_acceptable = sum(d["our_tm_score"] > 0.7 for d in data)
    ben_acceptable = sum(d["ben_best_tm_score"] > 0.7 for d in data)

    # Plot the high portion
    # ax.bar(x[0], our_high / total, label="High", color="blue")
    #
    # # Plot the acceptable portion ABOVE the high
    # ax.bar(x[0], (our_acceptable - our_high) / total, bottom=our_high / total, label="Acceptable", color="lightblue")

    our_high_pct = our_high / total
    our_acc_only_pct = (our_acceptable - our_high) / total

    ben_high_pct = ben_high / total
    ben_acc_only_pct = (ben_acceptable - ben_high) / total

    # our_high_pct = our_high / total * 100
    # ben_high_pct = ben_high / total * 100
    # our_lightblue = (our_acceptable - our_high) / total * 100
    # ben_lightblue = (ben_acceptable - ben_high) / total * 100


    # === PLOTTING ===
    # labels = ['Our TM-score', "Ben's TM-score"]
    # high = [our_high_pct, ben_high_pct]
    # light_blue = [our_lightblue, ben_lightblue]

    x = [0, 1]
    labels = ["Our models", "Ben's models"]

    fig, ax = plt.subplots()

    # High (blue)
    ax.bar(x[0], our_high_pct, color="blue")
    ax.bar(x[1], ben_high_pct, color="blue")

    # Acceptable only (light blue), stacked on top
    ax.bar(x[0], our_acc_only_pct, bottom=our_high_pct, color="lightblue", label="Acceptable (0.7-0.8)")
    ax.bar(x[1], ben_acc_only_pct, bottom=ben_high_pct, color="lightblue")

    # Annotations and legends
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of complexes")
    ax.set_title("TM-score comparison: Our models vs Ben's")
    ax.legend(["High (>0.8)", "Acceptable (0.7–0.8)"])

    plt.tight_layout()

    # === SAVE FIGURE ===
    output_dir = os.path.dirname(json_path)
    output_path = os.path.join(output_dir, "tm_score_quality_comparison.png")
    plt.savefig(output_path)
    print(f"✅ Plot saved to: {output_path}")



if __name__ == "__main__":
    # Run the check from the current directory
    #find_double_chain_proteins(base_dir="/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2")
    #main("/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2")
    plot_results_tm_score("/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2/tm_score_comparison.json")


