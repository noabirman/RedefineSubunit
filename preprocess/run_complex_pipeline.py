import os
import json
import sys
import subprocess
from string import ascii_uppercase


# --- 1. NEW FUNCTION TO CONVERT A FOLDER OF FASTAS TO ONE SUBUNIT JSON ---

def convert_folder_to_single_subunit_json(fasta_folder):
    """
    Converts all FASTA files in a folder into a single subunits_info.json format.

    Each FASTA file becomes a unique subunit, assigned a sequential alphabetical
    chain ID starting from 'A'.

    Args:
        fasta_folder (str): Path to the folder containing FASTA files.

    Returns:
        dict: The structured subunits JSON data for the entire complex.
    """
    complex_subunits = {}
    chain_generator = iter(ascii_uppercase)  # A, B, C, ...

    print("Collecting and validating FASTA files...")

    for filename in sorted(os.listdir(fasta_folder)):
        if filename.endswith(('.fasta', '.fa', '.fna')):
            fasta_path = os.path.join(fasta_folder, filename)

            try:
                # 1. FASTA Parsing
                with open(fasta_path, 'r') as f:
                    lines = f.readlines()

                header = ""
                sequence = ""

                # Simple parser: takes the first header and joins all sequence lines
                for line in lines:
                    line = line.strip()
                    if line.startswith('>'):
                        if not header:  # Only process the first header line
                            header = line.strip('>')
                    elif header:
                        sequence += line

                if not sequence:
                    print(f"  -> WARNING: Skipping {filename}. No sequence found.")
                    continue

                # 2. Assign IDs and Names
                # Subunit Key: Use the filename (without extension) as the unique subunit key
                subunit_key = filename.split('.')[0].upper()

                # Name: Use the full description from the FASTA header
                subunit_name = header.split('|')[-1].strip() if '|' in header else header

                # Chain Name: Get the next available single letter
                chain_id = next(chain_generator, None)
                if chain_id is None:
                    raise ValueError("Ran out of unique alphabetical chain IDs (A-Z).")

                # 3. Build Subunit Entry
                complex_subunits[subunit_key] = {
                    "name": subunit_name,
                    "chain_names": [chain_id],  # Assign a unique chain letter
                    "start_res": 1,
                    "sequence": sequence
                }
                print(f"  -> Added Subunit '{subunit_key}' (Chain '{chain_id}').")

            except Exception as e:
                print(f"  -> ERROR parsing {filename}: {e}. Skipping this file.")
                continue

    if not complex_subunits:
        raise RuntimeError("No valid FASTA files found or parsed in the folder.")

    return complex_subunits


# --- 2. MAIN EXECUTION LOGIC (REVISED) ---

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print("Usage: python3 run_complex_pipeline.py <folder_with_fastas>")
        sys.exit(1)

    fasta_folder = os.path.abspath(sys.argv[1])
    if not os.path.isdir(fasta_folder):
        print(f"Error: Folder '{fasta_folder}' not found.")
        sys.exit(1)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    preprocess_script_path = os.path.join(script_dir, "preprocess.py")

    # --- Define a SINGLE target output directory for the entire complex ---
    complex_dir_name = os.path.basename(fasta_folder).replace(" ", "_")
    output_dir = os.path.join(fasta_folder, f"AF3_COMPLEX_RESULTS")
    os.makedirs(output_dir, exist_ok=True)

    print(f"--- Starting Complex Preparation: {complex_dir_name} ---")

    # A. Convert ALL FASTAs in the folder to a SINGLE JSON
    try:
        subunit_json_data = convert_folder_to_single_subunit_json(fasta_folder)
        json_path = os.path.join(output_dir, "subunits_info.json")

        with open(json_path, 'w') as f:
            json.dump(subunit_json_data, f, indent=4)
        print(f"\n✅ Consolidated Subunit JSON created at: {json_path}")

    except Exception as e:
        print(f"\n❌ FATAL ERROR during JSON consolidation: {e}")
        sys.exit(1)

    # B. Run preprocess.py in DEFAULT Mode
    # NOTE: The script is called with only <dir_path> and <json_path>,
    # which triggers the "merge_duplicate_chain_sequences" and
    # "subunit_to_fasta" path (Default Mode).
    try:
        print("\n🚀 Running preprocess.py in Default (Subunit) Mode...")
        subprocess.run([
            "python3", preprocess_script_path,
            output_dir,  # <dir_path>
            json_path,  # <json_path>
        ], check=True)
        print(f"\n✨ Preprocessing complete for {complex_dir_name}. Check results in {output_dir}")

    except subprocess.CalledProcessError as e:
        print(f"\n❌ ERROR during preprocess.py execution. The process failed.")
        print(f"  Check the standard error output for details.")
        print(f"  Error: {e}")

    print("-" * 50)