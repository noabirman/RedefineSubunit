import os
import re
import sys
import json
import requests
import subprocess
from glob import glob
from typing import Tuple
import pandas as pd

MMALIGN_PATH = "/cs/labs/dina/bshor/scripts/MMalign"


def download_pdb_complex(pdb_id: str, out_dir: str = ".") -> str:
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    out_path = os.path.join(out_dir, f"{pdb_id}.pdb")
    if os.path.exists(out_path):
        return out_path
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(f"Failed to download PDB file: {pdb_id}")
    with open(out_path, "w") as f:
        f.write(response.text)
    return out_path


def get_tm_score_rmsd_mmalign(ref_path: str, sample_path: str) -> Tuple[float, float]:
    output = subprocess.check_output([MMALIGN_PATH, sample_path, ref_path]).decode()
    tm_scores = list(map(float, re.findall(r"TM-score= ([0-9]*[.]?[0-9]+)", output)))
    rmsds = list(map(float, re.findall(r"RMSD= *([0-9]*[.]?[0-9]+)", output)))
    assert len(tm_scores) == 2 and len(rmsds) >= 1
    return max(tm_scores), rmsds[0]



def process_complex(complex_dir: str, ben_scores: dict, variant_name: str):
    import os
    from glob import glob

    pdb_id = os.path.basename(complex_dir)
    print(f"\n🔍 Processing {pdb_id} [{variant_name}]")

    # Paths
    # if variant_name == "combfold_high":
    #     combfold_dir = os.path.join(complex_dir, "combfold", "results_high", "assembled_results")
    #     score_file = os.path.join(complex_dir, "combfold", "results_high", "tm_score.txt")
    # else:
    variant_dir = os.path.join(complex_dir, variant_name)
    combfold_dir = os.path.join(variant_dir, "results", "assembled_results")
    score_file = os.path.join(variant_dir, "tm_score.txt")

    # Ref structure
    ref_pdb = download_pdb_complex(pdb_id, complex_dir)

    if variant_name == "combfold_trivial":
        pattern = "output_clustered_0.pdb"
    else:
        pattern = "cb_*_output_0.pdb"

    clustered_models = sorted(glob(os.path.join(combfold_dir, pattern)))

    # Fallback for non-trivial variants if cb_* not found
    if not clustered_models and variant_name != "combfold_trivial":
        fallback_pattern = "output_clustered_0.pdb"
        clustered_models = sorted(glob(os.path.join(combfold_dir, fallback_pattern)))
        if clustered_models:
            print("✅ CombFold succeeded with output_clustered_0.pdb")

    if not clustered_models:
        print(f"⚠️  No clustered models found for {pdb_id} [{variant_name}]")
        return None

    # === CASE 1: If TM-score already computed, just read and return ===
    if os.path.exists(score_file):
        print(f"📄 Using existing TM-score file: {score_file}")
        best_model, best_tm, best_rmsd, ben_best_model, ben_best_tm = None, -1, -1, None, -1
        with open(score_file) as f:
            for line in f:
                if line.lower().startswith("our best model:"):
                    best_model = line.split(":")[1].strip()
                elif line.lower().startswith("tm-score:"):
                    best_tm = float(line.split(":")[1].strip())
                elif line.lower().startswith("rmsd:"):
                    best_rmsd = float(line.split(":")[1].strip())
                elif line.lower().startswith("ben's best model:"):
                    ben_best_model = line.split(":")[1].strip()
                elif "ben" in line.lower() and "tm-score" in line.lower():
                    ben_best_tm = float(line.split(":")[1].strip())
        return {
            "pdb_id": pdb_id.lower(),
            "scores": {
                variant_name: {"tm_score": best_tm, "rmsd": best_rmsd},
                "ben": {"tm_score": ben_best_tm} if ben_best_tm != -1 else {}
            }
        }

    # === CASE 2: Compute TM-score now ===

    # Align and compute TM-score
    best_tm, best_rmsd, best_model = -1, float("inf"), None
    for model in clustered_models:
        try:
            tm, rmsd = get_tm_score_rmsd_mmalign(ref_pdb, model)
            print(f"  → {os.path.basename(model)}: TM-score={tm:.5f}, RMSD={rmsd:.2f}")
            if tm > best_tm:
                best_tm = tm
                best_rmsd = rmsd
                best_model = os.path.basename(model)
        except Exception as e:
            print(f"  ❌ Error processing {model}: {e}")

    # Ben's scores
    pdb_id_upper = pdb_id.upper()
    ben_score_entries = ben_scores.get(pdb_id_upper, {}).get("scores", {})
    ben_best_tm = -1
    ben_best_model = None
    for model_name, score in ben_score_entries.items():
        tm_score = score.get("tm_score", -1)
        if tm_score > ben_best_tm:
            ben_best_tm = tm_score
            ben_best_model = model_name

    print(f"✅ Best of ours: {best_model} (TM={best_tm:.5f}, RMSD={best_rmsd:.2f})")
    if ben_best_model:
        print(f"✅ Best of Ben: {ben_best_model} (TM={ben_best_tm:.5f})")
    else:
        print("⚠️  No Ben results found for this PDB")

    # Write TM-score file
    with open(score_file, "w") as f:
        f.write(f"Our best model: {best_model}\n")
        f.write(f"TM-score: {best_tm:.5f}\n")
        f.write(f"RMSD: {best_rmsd:.2f}\n\n")
        if ben_best_model:
            f.write(f"Ben's best model: {ben_best_model}\n")
            f.write(f"TM-score: {ben_best_tm:.5f}\n")
        else:
            f.write("Ben's best model: Not found\n")

    # Return scores
    result = {
        "pdb_id": pdb_id.lower(),
        "scores": {
            variant_name: {"tm_score": best_tm, "rmsd": best_rmsd},
        }
    }
    if ben_best_tm != -1:
        result["scores"]["ben"] = {"tm_score": ben_best_tm}
    return result


def run_variant_on_all_complexes(root_dir: str, ben_scores: dict, variant_name: str, merged_results=None):
    if merged_results is None:
        merged_results = []

    merged_by_pdb = {entry["pdb_id"]: entry for entry in merged_results}

    for entry in sorted(os.listdir(root_dir)):
        complex_path = os.path.join(root_dir, entry)
        if os.path.isdir(complex_path) and os.path.isdir(os.path.join(complex_path, variant_name)):
            result = process_complex(complex_path, ben_scores, variant_name)
            if result:
                pdb_id = result["pdb_id"]
                if pdb_id not in merged_by_pdb:
                    merged_by_pdb[pdb_id] = {"pdb_id": pdb_id, "scores": {}}
                for k, v in result["scores"].items():
                    merged_by_pdb[pdb_id]["scores"][k] = v

    # Always include Ben’s scores if present
    for pdb_id, merged_entry in merged_by_pdb.items():
        if "ben" not in merged_entry["scores"]:
            ben_entry = ben_scores.get(pdb_id.upper(), {}).get("scores", {})
            if ben_entry:
                # Pick best Ben model TM-score
                best_ben_tm = max((v.get("tm_score", -1) for v in ben_entry.values()), default=-1)
                if best_ben_tm > -1:
                    merged_entry["scores"]["ben"] = {"tm_score": best_ben_tm}

    merged_results = list(merged_by_pdb.values())

    # Save to comparison JSON
    variant_out_path = os.path.join(root_dir, variant_name, "tm_score_comparison.json")
    os.makedirs(os.path.dirname(variant_out_path), exist_ok=True)
    with open(variant_out_path, "w") as f:
        json.dump(merged_results, f, indent=2)
    print(f"💾 Saved comparison results to {variant_out_path}")

    return merged_results

def plot_tm_score_summary(all_results, out_dir):
    import matplotlib.pyplot as plt
    import numpy as np
    #os.makedirs(out_dir, exist_ok=True)
    #output_path = os.path.join(root_dir, "all_tm_score_comparison.json")

    groups = list(all_results[0]["scores"].keys())
    tm_by_group = {g: [] for g in groups}

    for entry in all_results:
        for g in groups:
            tm_score = entry.get("scores", {}).get(g, {}).get("tm_score", None)
            if tm_score is not None:
                tm_by_group[g].append(tm_score)

    above_07 = [np.mean([s >= 0.7 for s in tm_by_group[g]]) * 100 for g in groups]
    above_08 = [np.mean([s >= 0.8 for s in tm_by_group[g]]) * 100 for g in groups]

    x = np.arange(len(groups))
    width = 0.35

    fig, ax = plt.subplots()
    ax.bar(x - width/2, above_07, width, label='TM ≥ 0.7')
    ax.bar(x + width/2, above_08, width, label='TM ≥ 0.8')

    ax.set_ylabel('Percentage of complexes')
    ax.set_title('TM-score thresholds by group')
    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=45, ha='right')
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "tm_score_summary.png"))
    plt.close()

def main_tm_score_script(root_dir: str, ben_json_path: str):
    with open(ben_json_path) as f:
        ben_scores = json.load(f)

    all_variants = [
        "combfold",
        "combfold_all",
        "combfold_trivial",
        "combfold_us_trivial",
        "combfold_high"
    ]

    all_results = {}
    for variant in all_variants:
        print(f"\n🧪 Running TM-score evaluation for: {variant}")
        all_results = run_variant_on_all_complexes(root_dir, ben_scores, variant, merged_results=all_results)
    # Save final comparison JSON
    output_path = os.path.join(root_dir, "tm_score_comparison.json")
    with open(output_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"💾 Final comparison JSON saved to: {output_path}")

    with open(output_path) as f:
        all_results = json.load(f)

    plots =os.path.join(root_dir,"plots")
    os.makedirs(plots, exist_ok=True)
    plot_tm_score_summary(all_results, out_dir=plots)


def check_complexes(root_dir):
    variants = [
        "combfold",
        "combfold_all",
        "combfold_trivial",
        "combfold_us_trivial",
        "combfold_high"
    ]

    records = []
    for complex_id in sorted(os.listdir(root_dir)):
        complex_dir = os.path.join(root_dir, complex_id)
        if not os.path.isdir(complex_dir):
            continue

        row = {"complex": complex_id}
        for variant in variants:
            variant_dir = os.path.join(complex_dir, variant, "results", "assembled_results")
            if os.path.exists(os.path.join(variant_dir, "output_clustered_0.pdb")):
                row[variant] = "success"
            else:
                row[variant] = "failed"
        records.append(row)

    df = pd.DataFrame(records)
    return df

if __name__ == "__main__":
    if len(sys.argv) != 3:
        # python3 tm_score_new.py /cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2 /cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2/combfold_results.json
        print("Usage: python run_all_tm_scores.py <DONE_MSA2_directory> <combfold_results.json>")
        sys.exit(1)

    root_dir = os.path.abspath(sys.argv[1])
    ben_json_path = os.path.abspath(sys.argv[2])
    main_tm_score_script(root_dir, ben_json_path)
    #check which of the complexs sucsess run fully
    df = check_complexes(root_dir)
    print(df.to_string(index=False))
    df.to_csv(os.path.join(root_dir, "clustered_presence_table.csv"), index=False)