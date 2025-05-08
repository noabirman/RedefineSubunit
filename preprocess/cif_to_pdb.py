import os
import sys
from Bio.PDB import MMCIFParser, PDBIO

def convert_cif_to_pdb(cif_path: str, pdb_path: str) -> None:
    """
    Convert a CIF file to PDB format.

    Parameters
    ----------
    cif_path : str
        Path to the input CIF file.
    pdb_path : str
        Path to the output PDB file.
    """
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(os.path.basename(cif_path), cif_path)
    except Exception as e:
        print(f"Error parsing CIF file '{cif_path}': {e}")
        return

    io = PDBIO()
    io.set_structure(structure)
    try:
        io.save(pdb_path)
        print(f"Successfully converted '{cif_path}' to '{pdb_path}'.")
    except Exception as e:
        print(f"Error writing PDB file '{pdb_path}': {e}")

def batch_convert_generic(input_dir: str, output_dir: str) -> None:
    """
    Convert all CIF files in the input directory to PDB format in the output directory.

    Parameters
    ----------
    input_dir : str
        Directory containing CIF files.
    output_dir : str
        Directory to save converted PDB files.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory '{output_dir}'.")

    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.cif'):
                cif_path = os.path.join(root, file)
                pdb_filename = os.path.splitext(file)[0] + '.pdb'
                pdb_path = os.path.join(output_dir, pdb_filename)
                convert_cif_to_pdb(cif_path, pdb_path)

def batch_convert_specific(input_dir: str, output_dir: str) -> None:
    """
    Convert specific CIF files within subdirectories of the input directory to PDB format.

    Parameters
    ----------
    input_dir : str
        Directory containing subdirectories with CIF files.
    output_dir : str
        Directory to save converted PDB files.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory '{output_dir}'.")

    for item in os.listdir(input_dir):
        item_path = os.path.join(input_dir, item)
        if os.path.isdir(item_path):
            cif_filename = f"{item}_model.cif"
            cif_path = os.path.join(item_path, cif_filename)

            if os.path.isfile(cif_path):
                pdb_filename = f"{item}.pdb"
                pdb_path = os.path.join(output_dir, pdb_filename)
                convert_cif_to_pdb(cif_path, pdb_path)
            else:
                print(f"Warning: CIF file '{cif_filename}' not found in '{item_path}'. Skipping.")

def main():
    # Verify we have enough arguments
    if len(sys.argv) < 3:
        print("Usage: python cif_to_pdb.py <input_dir> <output_dir> ")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    command = sys.argv[3] if len(sys.argv) > 3 else None
    batch_convert_generic(input_dir, output_dir)

if __name__ == '__main__':
    main()
