import os
import subprocess
if __name__ == '__main__':
    # List of complexes
    complexes = ['7pkn', '7zkq', '7oba', '8a5o', '7use', '7uic', '8cte']

    # Base directory containing all complex folders
    complex_dir = "/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2"

    # Optional: output folder name inside each complex (default in your script is 'models')
    output_folder_name = "models"

    for complex_name in complexes:
        parent_dir = os.path.join(complex_dir, complex_name, "combfold")

        # Make sure the complex folder exists
        if not os.path.exists(parent_dir):
            print(f"Complex folder not found: {parent_dir}, skipping.")
            continue

        # Call the existing cif_to_pdb.py script
        # Note: assumes cif_to_pdb.py is in the same folder as this script, otherwise provide full path
        subprocess.run([
            "python", "cif_to_pdb.py", parent_dir
        ])
