import os
def check_dup_and_size_missmatch(models_folder):
    if not os.path.exists(models_folder):
        print(f"  → Models folder does not exist, skipping: {models_folder}\n")
        return False  # No duplicates if folder doesn't exist
    # Path to the folder containing your PDB files
    #complex_folder = "/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2/7use/combfold/models"
    has_dup_with_diff_size = False
    # Collect all pdb filenames
    pdb_files = [f for f in os.listdir(models_folder) if f.endswith(".pdb")]

    # Dictionary to store normalized pair -> filename
    pairs = {}

    for f in pdb_files:
        # Example: "a_b_model.pdb"
        if not f.endswith("_model.pdb"):
            continue

        base = f[:-10]  # remove "_model.pdb"
        parts = base.split("_")
        if len(parts) != 2:
            continue  # skip unexpected file name patterns

        # Normalize so a_b and b_a become the same key
        pair_key = tuple(sorted(parts))

        if pair_key in pairs:
            pairs[pair_key].append(f)
        else:
            pairs[pair_key] = [f]

    # Now check for duplicates and compare file sizes
    dup_pairs = 0
    for pair, files in pairs.items():
        if len(files) > 1:
            print(f"Found duplicate pair: {pair}")
            dup_pairs = dup_pairs + 1
            for fname in files:
                path = os.path.join(models_folder, fname)
                size = os.path.getsize(path)
                print(f"  {fname}: {size} bytes")

            sizes = [os.path.getsize(os.path.join(models_folder, f)) for f in files]
            if len(set(sizes)) == 1:
                print("  → Same file size\n")
            else:
                print("  → Different file sizes\n")
                has_dup_with_diff_size = True
    print(f"Total number of duplicate pairs: {dup_pairs}")
    return has_dup_with_diff_size

if __name__ == "__main__":
    COMPLEXES = ["8a3t", "8adl", "8cte", "8f5o", "7wff", "7e8t", "8hil", "7t3b", "7oba", "7uic", "7pkn", "7xho", "7zkq",
                 "8a5o", "8fee", "8bel", "7qve", "7arc", "7ozn", "8adn", "7t2r", "7p2y", "7qru", "7use", "8e9g"]
    complexes_root = "/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2"
    complex_folders = [os.path.join(complexes_root, d) for d in os.listdir(complexes_root)
                       if os.path.isdir(os.path.join(complexes_root, d))]

    problematic_folders = []
    for complex_folder in complex_folders:
        complex_name = os.path.basename(complex_folder)
        if(complex_name in COMPLEXES):
            print(f"checking {complex_name}....")
            models_folder = os.path.join(complex_folder, "combfold", "models")
            if(check_dup_and_size_missmatch(models_folder)):
                problematic_folders.append(complex_name)
    print(problematic_folders)