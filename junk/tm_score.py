import gzip
import re
import shutil
import subprocess
import sys
import argparse
from Bio.PDB import PDBList
import os

# def download_pdb_complex(pdb_id: str, out_dir: str = ".") -> str:
#     """
#     Downloads the full PDB structure (complex) in PDB format.
#
#     Parameters:
#         pdb_id (str): PDB identifier (e.g., "1A4U").
#         out_dir (str): Directory to save the downloaded file.
#
#     Returns:
#         str: Path to the downloaded PDB file.
#     """
#     pdb_id = pdb_id.lower()
#     pdbl = PDBList(verbose=0)
#     pdb_gz_path = pdbl.retrieve_pdb_file(pdb_code=pdb_id, file_format="pdb", pdir=out_dir)
#
#     # Unzip if the file is gzipped (default behavior is to download .ent.gz)
#     if pdb_gz_path.endswith(".gz"):
#         import gzip, shutil
#         pdb_path = pdb_gz_path[:-3]  # remove .gz
#         with gzip.open(pdb_gz_path, 'rb') as f_in, open(pdb_path, 'wb') as f_out:
#             shutil.copyfileobj(f_in, f_out)
#         os.remove(pdb_gz_path)  # optionally delete the .gz file
#     else:
#         pdb_path = pdb_gz_path
#
#     return pdb_path
# def download_pdb_complex(pdb_id: str, out_dir: str = ".") -> str:
#     """
#     Downloads the full PDB structure (complex) in PDB format and renames to <pdb_id>.pdb.
#
#     Parameters:
#         pdb_id (str): PDB identifier (e.g., "7T3B").
#         out_dir (str): Directory to save the downloaded file.
#
#     Returns:
#         str: Path to the downloaded and renamed PDB file.
#     """
#     pdb_id = pdb_id.lower()
#     pdbl = PDBList(verbose=0)
#
#     # Force download into flat directory
#     pdb_gz_path = pdbl.retrieve_pdb_file(pdb_code=pdb_id, file_format="pdb", pdir=out_dir, overwrite=True)
#
#     # Biopython saves as something like: <out_dir>/pdb/pdbXXXX.ent.gz
#     # Uncompress and rename to XXXX.pdb
#     if pdb_gz_path.endswith(".gz"):
#         pdb_filename = f"{pdb_id}.pdb"
#         pdb_path = os.path.join(out_dir, pdb_filename)
#
#         with gzip.open(pdb_gz_path, "rb") as f_in, open(pdb_path, "wb") as f_out:
#             shutil.copyfileobj(f_in, f_out)
#
#         os.remove(pdb_gz_path)  # Clean up .gz
#         return pdb_path
#
#     else:
#         return pdb_gz_path
import requests

# def download_pdb_complex(pdb_id: str, out_dir: str = ".") -> str:
#     """
#     Downloads the PDB file (not mmCIF) for a given PDB ID using RCSB HTTP API.
#
#     Parameters:
#         pdb_id (str): PDB identifier (e.g., "7T3B")
#         out_dir (str): Directory to save the file
#
#     Returns:
#         str: Path to the downloaded .pdb file
#     """
#     pdb_id = pdb_id.lower()
#     url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
#     out_path = os.path.join(out_dir, f"{pdb_id}.pdb")
#
#     print(f"Downloading PDB from: {url}")
#     response = requests.get(url)
#     if response.status_code != 200:
#         raise Exception(f"Failed to download PDB file: {pdb_id}")
#
#     with open(out_path, "w") as f:
#         f.write(response.text)
#
#     print(f"Saved PDB to: {out_path}")
#     return out_path
#
#
# def get_tm_score_rmsd_mmalign(ref_complex_path: str, sample_complex_path: str):
#   MMALIGN_PATH = "/cs/labs/dina/bshor/scripts/MMalign"
#   mm_output = subprocess.check_output([MMALIGN_PATH, sample_complex_path, ref_complex_path]).decode()
#   tm_scores = list(map(float, re.findall(r"TM-score= ([0-9]*[.]?[0-9]+)", mm_output)))
#   rmsd = list(map(float, re.findall(r"RMSD= *([0-9]*[.]?[0-9]+)", mm_output)))
#   assert len(tm_scores) == 2 and len(rmsd) == 1
#   return max(tm_scores), rmsd[0]
#
#
# def tm_score_for_complex():
#     parser = argparse.ArgumentParser(description="Process CombFold result folder.")
#     parser.add_argument("folder", help="Path to the folder named after the PDB ID")
#
#     args = parser.parse_args()
#     given_folder = os.path.abspath(args.folder)
#     pdb_id = os.path.basename(given_folder)
#     output_txt_path = os.path.join(given_folder, "tm_score.txt")
#
#
#     # Construct path to output_clustered_0.pdb
#     model_path = os.path.join(
#         given_folder,
#         "combfold",
#         "results",
#         "assembled_results",
#         "output_clustered_0.pdb"
#     )
#     # Construct full path to reference PDB
#     ref_pdb_path = os.path.abspath(os.path.join(given_folder, pdb_id + ".pdb"))
#     print(f"PDB ID: {pdb_id}")
#     print(f"Expected output path: {given_folder}")
#
#     # Download PDB reference (entire complex in PDB format)
#     download_pdb_complex(pdb_id, given_folder)
#     tm_score, rmsd = get_tm_score_rmsd_mmalign(ref_pdb_path, model_path)
#     print(f"TM-score: {tm_score}, RMSD: {rmsd}")
#     with open(output_txt_path, "w") as f:
#         f.write(f"TM-score: {tm_score}\n")
#         f.write(f"RMSD: {rmsd}\n")
#     print(f"Saved scores to {output_txt_path}")
#
import os
import re
import json
import argparse
import subprocess
import requests
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


def process_complex(complex_dir: str, ben_scores: dict):
    pdb_id = os.path.basename(complex_dir)
    print(f"\n🔍 Processing {pdb_id}")

    combfold_dir = os.path.join(complex_dir, "combfold", "results", "assembled_results")
    ref_pdb = download_pdb_complex(pdb_id, complex_dir)

    clustered_models = sorted(glob(os.path.join(combfold_dir, "cb_*_output_0.pdb")))
    if not clustered_models:
        print(f"⚠️  No clustered models found for {pdb_id}")
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
    pdb_id = os.path.basename(complex_dir).upper()
    ben_score_entries = ben_scores.get(pdb_id, {}).get("scores", {})
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
    score_file = os.path.join(complex_dir, "tm_score.txt")
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


def main(root_dir: str, ben_json_path: str):
    with open(ben_json_path) as f:
        ben_scores = json.load(f)

    results = []
    for entry in sorted(os.listdir(root_dir)):
        complex_path = os.path.join(root_dir, entry)
        if os.path.isdir(complex_path):
            result = process_complex(complex_path, ben_scores)
            if result:
                results.append(result)

    out_path = os.path.join(root_dir, "tm_score_comparison.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved comparison results to {out_path}")


import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_all_tm_scores.py <DONE_MSA2_directory> <combfold_results.json>")
        sys.exit(1)

    root_dir = os.path.abspath(sys.argv[1])
    ben_json_path = os.path.abspath(sys.argv[2])
    main(root_dir, ben_json_path)
#
# if __name__ == "__main__":
#     tm_score_for_complex("/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2/7arc")





