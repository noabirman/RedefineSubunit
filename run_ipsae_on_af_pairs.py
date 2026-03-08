import os
import glob
import subprocess
import shutil
import sys

# Path to the ipsae.py tool
IPSAE_SCRIPT = "/cs/labs/dina/noabirman/RedefineSubunit/ipsae.py"


def run_pipeline(root_dir, pae_cutoff, dist_cutoff):

    # Create output directory name dynamically
    output_dir = os.path.join(
    os.path.dirname(root_dir),
    f"ipsae_results_{pae_cutoff:02d}_{dist_cutoff:02d}"
)
    print("Will create output folder:", output_dir)

    # 1. Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # 2. Iterate through all subdirectories in the root folder
    for folder_name in os.listdir(root_dir):
        folder_path = os.path.join(root_dir, folder_name)

        # Skip files, process only directories
        if not os.path.isdir(folder_path):
            continue

        print(f"Processing folder: {folder_name}...")

        # 3. Find required files
        cif_files = glob.glob(os.path.join(folder_path, "*.cif"))
        json_candidates = glob.glob(os.path.join(folder_path, "*confidences.json"))

        full_json_files = [
            f for f in json_candidates
            if "summary" not in os.path.basename(f)
        ]

        if not cif_files or not full_json_files:
            print(f"  [Skipping] Missing .cif or confidences.json in {folder_name}")
            continue

        structure_file = cif_files[0]
        pae_file = full_json_files[0]

        # 4. Construct command
        cmd = [
            sys.executable, IPSAE_SCRIPT,
            pae_file,
            structure_file,
            str(pae_cutoff),
            str(dist_cutoff)
        ]

        try:
            subprocess.run(cmd, check=True)
            print(f"  -> Successfully ran ipsae for {folder_name}")

            # 5. Move generated files
            cif_filename = os.path.basename(structure_file)
            stem = cif_filename.replace(".cif", "")

            # pae_str = str(pae_cutoff)
            # dist_str = str(dist_cutoff)
            pae_str = f"{pae_cutoff:02d}"
            dist_str = f"{dist_cutoff:02d}"

            generated_stem = f"{stem}_{pae_str}_{dist_str}"
            extensions = [".txt", "_byres.txt", ".pml"]

            for ext in extensions:
                filename = generated_stem + ext
                src_path = os.path.join(folder_path, filename)
                dst_path = os.path.join(
                    output_dir,
                    f"{folder_name}_{filename}"
                )

                if os.path.exists(src_path):
                    shutil.move(src_path, dst_path)

            print("  -> Moved outputs to results folder.")

        except subprocess.CalledProcessError as e:
            print(f"  [Error] Failed to run ipsae for {folder_name}: {e}")
        except Exception as e:
            print(f"  [Error] An issue occurred: {e}")


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage:")
        print("python3 run_ipsae_on_af_pairs.py <ROOT_DIR> <PAE_CUTOFF> <DIST_CUTOFF>")
        sys.exit(1)

    root_dir = sys.argv[1]
    pae_cutoff = int(sys.argv[2])
    dist_cutoff = int(sys.argv[3])

    run_pipeline(root_dir, pae_cutoff, dist_cutoff)