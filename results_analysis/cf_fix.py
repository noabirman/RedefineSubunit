import subprocess
import os
if __name__ == '__main__':
    # List of complexes
    complexes = ['7pkn', '7zkq', '7oba', '8a5o', '7use', '7uic', '8cte']

    # Base directory containing all complex folders
    base_dir = "/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2"

    # Path to your run_cf.sh script
    run_script = "/cs/labs/dina/noabirman/RedefineSubunit/combfold/run_cf.sh"


    for complex_name in complexes:
        complex_path = os.path.join(base_dir, complex_name)

        if not os.path.exists(complex_path):
            print(f"Complex folder not found: {complex_path}, skipping.")
            continue

        # Run the command: ./run_cf.sh <complex_path> 1
        subprocess.run([run_script, complex_path, "1"])
