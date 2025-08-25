import os
import sys
from Bio.PDB import MMCIFParser, PDBIO

def convert_cif_to_pdb(cif_path: str, pdb_path: str) -> None:
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
        # print(f"Successfully converted '{os.path.basename(cif_path)}' to '{os.path.basename(pdb_path)}'.")
    except Exception as e:
        print(f"Error writing PDB file '{pdb_path}': {e}")

def batch_convert_generic(input_dir: str, output_dir: str) -> None:

    files_converted = 0

    for subdir in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, subdir)
        if os.path.isdir(subdir_path):
            for file in os.listdir(subdir_path):
                if file.endswith('.cif'):
                    cif_path = os.path.join(subdir_path, file)
                    pdb_filename = os.path.splitext(file)[0] + '.pdb'
                    pdb_path = os.path.join(output_dir, pdb_filename)
                    convert_cif_to_pdb(cif_path, pdb_path)
                    files_converted += 1

    print(f"Done converting {files_converted} CIF files.")

if __name__ == '__main__':
    if len(sys.argv) not in [2, 3]:
        print("Usage: python cif_to_pdb.py <input_parent_dir> [output_dir]")
        sys.exit(1)

    parent_dir = os.path.abspath(sys.argv[1])
    input_dir = os.path.join(parent_dir, 'af_pairs')
    if len(sys.argv) == 3:
        output_dir = os.path.abspath(sys.argv[2])
    else:
        output_dir = os.path.join(parent_dir, 'models')

    print("Converting CIF files to PDB files...")
    print(f"Input directory: {input_dir}")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: '{output_dir}'.")
    else:
        print(f"Output directory: {output_dir}")

    batch_convert_generic(input_dir, output_dir)