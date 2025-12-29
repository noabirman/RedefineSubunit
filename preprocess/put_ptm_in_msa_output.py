import json
import os


def process_msa_file(msa_json_path, phospho_sites, output_dir):
    """
    Process one A_D JSON file: adds phospho modifications and writes to output_dir.

    Args:
        msa_json_path (str): Path to the input JSON file.
        phospho_sites (dict): Phospho sites dictionary.
        output_dir (str): Path to output folder.
    """
    input_base = os.path.basename(msa_json_path)
    name, ext = os.path.splitext(input_base)
    if name.endswith("_data"):
        name = name[:-5]
    output_name = f"{name.upper()}{ext}"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_name)

    # load the JSON
    with open(msa_json_path, "r") as f:
        msa_file = json.load(f)

    # add phospho modifications
    for seq_entry in msa_file["sequences"]:
        protein = seq_entry["protein"]
        pid = protein["id"]
        ptm_list = [
            {"ptmType": p["ccd"], "ptmPosition": p["position"]}
            for p in phospho_sites.get(pid, [])
        ]
        protein["modifications"] = ptm_list

    # write output file
    with open(output_path, "w") as f:
        json.dump(msa_file, f, indent=2)

    print(f"Processed {msa_json_path} -> {output_path}")


def main(msa_output_root, phospho_file):
    # load phospho sites once
    with open(phospho_file, "r") as f:
        phospho_sites = json.load(f)

    # iterate over all subdirectories in msa_output_root
    for subdir in os.listdir(msa_output_root):
        subdir_path = os.path.join(msa_output_root, subdir)
        if not os.path.isdir(subdir_path):
            continue

        # find JSON file ending with "_data.json"
        json_files = [f for f in os.listdir(subdir_path) if f.endswith("_data.json")]
        if not json_files:
            print(f"No '_data.json' file found in {subdir_path}, skipping...")
            continue

        msa_json_path = os.path.join(subdir_path, json_files[0])

        # determine output directory (same logic as before)
        output_dir = os.path.join(os.path.dirname(msa_output_root), "af_triv_input")

        # process the JSON file
        process_msa_file(msa_json_path, phospho_sites, output_dir)

def run_on_one_file(phospho_file):
    with open(PHOSPHO_FILE, "r") as f:
        phospho_sites = json.load(f)
    msa_json_path ="/cs/labs/dina/noabirman/tcellsUniprots/AF3_COMPLEX_RESULTS/combfold/msa_output/a_ab/a_ab_data.json"
    output_dir = "/cs/labs/dina/noabirman/tcellsUniprots/AF3_COMPLEX_RESULTS/combfold/msa_output/a_ab/a_ab_try_pho.json"
    process_msa_file(msa_json_path, phospho_sites, output_dir)

if __name__ == "__main__":
    MSA_OUTPUT_ROOT = "/cs/labs/dina/noabirman/tcellsUniprots/AF3_COMPLEX_RESULTS/combfold/msa_output" #before was "/cs/labs/dina/noabirman/tcellsUniprots/AF3_COMPLEX_RESULTS/msa_output"
    PHOSPHO_FILE = "phospho_final.json" #before was phospho_after_preprocess.json
    main(MSA_OUTPUT_ROOT, PHOSPHO_FILE)

    #run_on_one_file(PHOSPHO_FILE)
