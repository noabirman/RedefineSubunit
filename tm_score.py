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

def download_pdb_complex(pdb_id: str, out_dir: str = ".") -> str:
    """
    Downloads the PDB file (not mmCIF) for a given PDB ID using RCSB HTTP API.

    Parameters:
        pdb_id (str): PDB identifier (e.g., "7T3B")
        out_dir (str): Directory to save the file

    Returns:
        str: Path to the downloaded .pdb file
    """
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    out_path = os.path.join(out_dir, f"{pdb_id}.pdb")

    print(f"Downloading PDB from: {url}")
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(f"Failed to download PDB file: {pdb_id}")

    with open(out_path, "w") as f:
        f.write(response.text)

    print(f"Saved PDB to: {out_path}")
    return out_path


def get_tm_score_rmsd_mmalign(ref_complex_path: str, sample_complex_path: str):
  MMALIGN_PATH = "/cs/labs/dina/bshor/scripts/MMalign"
  mm_output = subprocess.check_output([MMALIGN_PATH, sample_complex_path, ref_complex_path]).decode()
  tm_scores = list(map(float, re.findall(r"TM-score= ([0-9]*[.]?[0-9]+)", mm_output)))
  rmsd = list(map(float, re.findall(r"RMSD= *([0-9]*[.]?[0-9]+)", mm_output)))
  assert len(tm_scores) == 2 and len(rmsd) == 1
  return max(tm_scores), rmsd[0]


def main():
    parser = argparse.ArgumentParser(description="Process CombFold result folder.")
    parser.add_argument("folder", help="Path to the folder named after the PDB ID")

    args = parser.parse_args()
    given_folder = os.path.abspath(args.folder)
    pdb_id = os.path.basename(given_folder)

    # Construct path to output_clustered_0.pdb
    model_path = os.path.join(
        given_folder,
        "combfold",
        "results",
        "assembled_results",
        "output_clustered_0.pdb"
    )
    # Construct full path to reference PDB
    ref_pdb_path = os.path.abspath(os.path.join(given_folder, pdb_id + ".pdb"))
    print(f"PDB ID: {pdb_id}")
    print(f"Expected output path: {given_folder}")

    # Download PDB reference (entire complex in PDB format)
    download_pdb_complex(pdb_id, given_folder)
    tm_score, rmsd = get_tm_score_rmsd_mmalign(ref_pdb_path, model_path)
    print(f"TM-score: {tm_score}, RMSD: {rmsd}")


if __name__ == "__main__":
    main()
