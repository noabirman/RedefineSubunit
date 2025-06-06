import os
import re
import sys
import json
import requests
import subprocess
from glob import glob
from typing import Tuple

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
    pdb_id = os.path.basename(complex_dir)
    print(f"\n🔍 Processing {pdb_id} [{variant_name}]")

    variant_dir = os.path.join(complex_dir, variant_name)
    combfold_dir = os.path.join(variant_dir, "results", "assembled_results")
    ref_pdb = download_pdb_complex(pdb_id, complex_dir)

    if "trivial" in variant_name.lower():
        pattern = "output_clustered_0.pdb"
    else:
        pattern = "cb_*_output_0.pdb"

    clustered_models = sorted(glob(os.path.join(combfold_dir, pattern)))

    if not clustered_models:
        print(f"⚠️  No clustered models found for {pdb_id} [{variant_name}]")
        return None

    best_tm, best_rmsd, best_model = -1, float("inf"), None
    for model in clustered_models:
        try:
            tm, rmsd = get_tm_score_rmsd_mmalign(ref_pdb, model)
            print(f"  → {os.path.basename(model)}: TM-score={tm:.5f}, RMSD={rmsd:.2f}")
            if tm > best_tm:
                best_tm, best_rmsd, best_model = tm, rmsd, os.path.basename(model)
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

    # Save to tm_score.txt
    score_file = os.path.join(variant_dir, "tm_score.txt")
    with open(score_file, "w") as f:
        f.write(f"Our best model: {best_model}\n")
        f.write(f"TM-score: {best_tm:.5f}\n")
        f.write(f"RMSD: {best_rmsd:.2f}\n\n")
        if ben_best_model:
            f.write(f"Ben's best model: {ben_best_model}\n")
            f.write(f"TM-score: {ben_best_tm:.5f}\n")
        else:
            f.write("Ben's best model: Not found\n")

    return {
        "pdb_id": pdb_id.lower(),
        "best_model": best_model,
        "our_tm_score": best_tm,
        "our_rmsd": best_rmsd,
        "ben_best_model": ben_best_model,
        "ben_best_tm_score": ben_best_tm
    }


def run_variant_on_all_complexes(root_dir: str, ben_scores: dict, variant_name: str):
    results = []
    for entry in sorted(os.listdir(root_dir)):
        complex_path = os.path.join(root_dir, entry)
        if os.path.isdir(complex_path) and os.path.isdir(os.path.join(complex_path, variant_name)):
            result = process_complex(complex_path, ben_scores, variant_name)
            if result:
                results.append(result)

    # Save to tm_score_comparison.json in the variant directory
    variant_out_path = os.path.join(root_dir, variant_name, "tm_score_comparison.json")
    os.makedirs(os.path.dirname(variant_out_path), exist_ok=True)
    with open(variant_out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"💾 Saved comparison results to {variant_out_path}")


def main(root_dir: str, ben_json_path: str):
    with open(ben_json_path) as f:
        ben_scores = json.load(f)

    variant_names = ["combfold_all", "combfold_trivial", "combfold_us_trivial"]
    for variant in variant_names:
        print(f"\n🧪 Running TM-score evaluation for: {variant}")
        run_variant_on_all_complexes(root_dir, ben_scores, variant)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        # python3 tm_score_new.py /cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2 /cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2/combfold_results.json
        print("Usage: python run_all_tm_scores.py <DONE_MSA2_directory> <combfold_results.json>")
        sys.exit(1)

    root_dir = os.path.abspath(sys.argv[1])
    ben_json_path = os.path.abspath(sys.argv[2])
    main(root_dir, ben_json_path)