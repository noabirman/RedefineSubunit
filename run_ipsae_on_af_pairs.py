import os
import glob
import subprocess
import shutil
import sys

# --- CONFIGURATION ---
# The folder containing your pair subfolders (e.g., current directory)
ROOT_DIR = "/cs/labs/dina/noabirman/tcellsUniprots/AF3_COMPLEX_RESULTS/af_pairs"

# The folder where you want to gather all outputs (created if not exists)
OUTPUT_DIR = "/cs/labs/dina/noabirman/tcellsUniprots/AF3_COMPLEX_RESULTS/ipsae_results"

# Cutoff parameters
PAE_CUTOFF = 15
DIST_CUTOFF = 10

# Path to the ipsae.py tool
IPSAE_SCRIPT = "/cs/labs/dina/noabirman/RedefineSubunit/ipsae.py"


# ---------------------

def run_pipeline():
    # 1. Create the output directory if it doesn't exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created output directory: {OUTPUT_DIR}")

    # 2. Iterate through all subdirectories in the root folder
    for folder_name in os.listdir(ROOT_DIR):
        folder_path = os.path.join(ROOT_DIR, folder_name)

        # Skip files, process only directories (and skip the output dir itself)
        if not os.path.isdir(folder_path) or folder_name == OUTPUT_DIR:
            continue

        print(f"Processing folder: {folder_name}...")

        # 3. Find the required files (.cif and _confidences.json)
        # We look for any .cif file and any JSON ending in 'confidences.json'
        cif_files = glob.glob(os.path.join(folder_path, "*.cif"))
        json_candidates = glob.glob(os.path.join(folder_path, "*confidences.json"))

        full_json_files = [f for f in json_candidates if "summary" not in os.path.basename(f)]

        if not cif_files or not full_json_files:
            print(f"  [Skipping] Missing .cif or confidences.json in {folder_name}")
            continue

        # Take the first match found
        structure_file = cif_files[0]
        pae_file = full_json_files[0]

        # 4. Construct the command
        # Syntax: python ipsae.py <pae_file> <structure_file> <pae_cutoff> <dist_cutoff>
        cmd = [
            sys.executable, IPSAE_SCRIPT,
            pae_file,
            structure_file,
            str(PAE_CUTOFF),
            str(DIST_CUTOFF)
        ]

        try:
            # 5. Run ipsae.py
            subprocess.run(cmd, check=True)
            print(f"  -> Successfully ran ipsae for {folder_name}")

            # 6. Move the output files
            # ipsae.py names outputs based on the input structure filename
            # Structure: name.cif -> Outputs: name_10_10.txt, name_10_10_byres.txt, etc.

            # Identify the stem used by ipsae logic
            cif_filename = os.path.basename(structure_file)
            stem = cif_filename.replace(".cif", "")

            # Format cutoffs as 2 digits (e.g., 5 -> "05", 10 -> "10")
            pae_str = f"{PAE_CUTOFF:02d}"
            dist_str = f"{DIST_CUTOFF:02d}"

            generated_stem = f"{stem}_{pae_str}_{dist_str}"

            # List of extensions ipsae generates
            extensions = [".txt", "_byres.txt", ".pml"]

            for ext in extensions:
                filename = generated_stem + ext
                src_path = os.path.join(folder_path, filename)
                dst_path = os.path.join(OUTPUT_DIR,
                                        f"{folder_name}_{filename}")  # Prepend folder name to avoid collisions

                if os.path.exists(src_path):
                    shutil.move(src_path, dst_path)

            print("  -> Moved outputs to results folder.")

        except subprocess.CalledProcessError as e:
            print(f"  [Error] Failed to run ipsae for {folder_name}: {e}")
        except Exception as e:
            print(f"  [Error] An issue occurred: {e}")


if __name__ == "__main__":
    run_pipeline()