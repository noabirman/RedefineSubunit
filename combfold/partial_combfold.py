COMBFOLD_PATH = '/cs/labs/dina/tsori/af3_example'
import os
import sys
sys.path.insert(0, COMBFOLD_PATH)
from CombFold.scripts.libs.prepare_complex import create_complexes

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <results_path>")
        sys.exit(1)

    results_path = os.path.abspath(sys.argv[1])
    output_folder = os.path.join(results_path, "assembled_results")
    assembly_dir = os.path.join(results_path, "_unified_representation", "assembly_output")

    print(f"Starting assembly process with results path: {results_path}")

    # Find the clusters file
    clusters_path = None
    for fname in os.listdir(assembly_dir):
        if fname.startswith("cb") and fname.endswith(".res"):
            clusters_path = os.path.join(assembly_dir, fname)
            break

    if clusters_path is None:
        print(f"Could not find clusters file (cb*.res) in {assembly_dir}, exiting")
        sys.exit(1)

    # Example usage of create_complexes (adjust arguments as needed)
    assembled_file = create_complexes(
        clusters_path,
        first_result=0,
        last_result=10,  # or another value
        output_folder=output_folder,
        output_cif=False  # or False, as needed
    )
    print(f"Assembly completed. Results saved to: {assembled_file}")
