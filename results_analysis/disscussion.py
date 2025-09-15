
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder
# import os
# import json


# from Bio.PDB import PDBParser, MMCIFParser, PPBuilder
# from Bio import pairwise2
# import os
# import json
#
# def extract_chain_sequence(structure_file, chain_id):
#     """Extract observed amino acid sequence for a chain from PDB or CIF file."""
#     ext = os.path.splitext(structure_file)[1].lower()
#
#     if ext in [".pdb", ".ent"]:
#         parser = PDBParser(QUIET=True)
#     elif ext in [".cif", ".mmcif"]:
#         parser = MMCIFParser(QUIET=True)
#     else:
#         raise ValueError(f"Unsupported file format: {ext}")
#
#     structure = parser.get_structure("structure", structure_file)
#     model = structure[0]
#     if chain_id not in model:
#         raise ValueError(f"Chain {chain_id} not found in {structure_file}")
#
#     chain = model[chain_id]
#     ppb = PPBuilder()
#     peptides = ppb.build_peptides(chain)
#
#     if not peptides:
#         return ""  # empty chain
#
#     return "".join([str(p.get_sequence()) for p in peptides])
#
#
# def reconstruct_expected_chain_sequence(subunit_info, chain_id):
#     """Reconstructs the full expected sequence for a chain from multiple subunits."""
#     fragments = []
#     for name, info in subunit_info.items():
#         if chain_id in info["chain_names"]:
#             start = info.get("start_res", 1)
#             fragments.append((start, info["sequence"], name))
#
#     if not fragments:
#         raise ValueError(f"No subunits found for chain {chain_id}")
#
#     # Sort by residue start
#     fragments.sort(key=lambda x: x[0])
#
#     expected_seq = ""
#     mapping = []  # (res_index, subunit_name)
#     for start, seq, subunit in fragments:
#         for i, aa in enumerate(seq):
#             res_index = start + i
#             expected_seq += aa
#             mapping.append((res_index, subunit))
#
#     return expected_seq, mapping
#
#

def extract_chain_sequence(structure_file, chain_id):
    """Extract observed amino acid sequence for a chain from PDB or CIF file."""
    ext = os.path.splitext(structure_file)[1].lower()

    if ext in [".pdb", ".ent"]:
        parser = PDBParser(QUIET=True)
    elif ext in [".cif", ".mmcif"]:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported file format: {ext}")

    structure = parser.get_structure("structure", structure_file)
    model = structure[0]
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found in {structure_file}")

    chain = model[chain_id]
    ppb = PPBuilder()
    peptides = ppb.build_peptides(chain)

    if not peptides:
        return ""  # empty chain

    return "".join([str(p.get_sequence()) for p in peptides])
def reconstruct_expected_chain_sequence(subunit_info, chain_id):
    """Reconstructs the full expected sequence for a chain from multiple subunits."""
    parts = []
    for name, info in subunit_info.items():
        if chain_id in info["chain_names"]:
            start = info.get("start_res", 1)
            seq = info["sequence"]
            parts.append((start, seq, name))

    if not parts:
        raise ValueError(f"No subunits found for chain {chain_id}")

    # Sort fragments by starting residue
    parts.sort(key=lambda x: x[0])

    # Build full sequence and map residues back to subunit origin
    expected_seq = ""
    mapping = []  # (res_index, subunit_name)
    for start, seq, subunit_name in parts:
        for i, aa in enumerate(seq):
            expected_seq += aa
            mapping.append((start + i, subunit_name))

    return expected_seq, mapping


import os
from Bio.PDB import PDBParser, is_aa
import json
import numpy as np
import pandas as pd

# Custom 3-letter → 1-letter map
aa_three_to_one = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O"
}

def get_pdb_residues(pdb_file, chain_id):
    """Extract residue numbers and 1-letter codes from a chain in the PDB."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = structure[0]
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found in PDB file")
    chain = model[chain_id]

    residues = {}
    for res in chain:
        if is_aa(res, standard=True):
            resseq = res.id[1]
            resname = res.resname
            aa = aa_three_to_one.get(resname, "X")
            residues[resseq] = aa
    return residues


def check_missing_and_mismatches(pdb_file, subunit_info, chain_id):
    """Check missing and mismatched residues for a given chain across all subunits."""
    pdb_residues = get_pdb_residues(pdb_file, chain_id)
    report = {}

    for sub_name, info in subunit_info.items():
        if chain_id not in info["chain_names"]:
            continue

        seq = info["sequence"]
        start = info.get("start_res", 1)
        expected = {start + i: aa for i, aa in enumerate(seq)}

        missing = []
        mismatches = []
        for resnum, aa in expected.items():
            if resnum not in pdb_residues:
                missing.append(resnum)
            elif pdb_residues[resnum] != aa:
                mismatches.append((resnum, aa, pdb_residues[resnum]))

        report[sub_name] = {
            "expected_range": f"{start}-{start+len(seq)-1}",
            "missing": missing,
            "mismatches": mismatches,
        }

    return report


def analyze_all_chains(pdb_file, subunit_info):
    """Run the check for each chain ID that appears in the JSON."""
    # Collect all chain IDs, including duplicates across subunits
    all_chains = set()
    for info in subunit_info.values():
        for chain_id in info["chain_names"]:
            all_chains.add(chain_id)

    results = {}
    for chain_id in sorted(all_chains):
        try:
            results[chain_id] = check_missing_and_mismatches(pdb_file, subunit_info, chain_id)
        except ValueError as e:
            results[chain_id] = {"error": str(e)}

    return results

import pandas as pd  # for table creation

def build_summary_table(results, subunit_info):
    rows = []

    # Gather all chain IDs
    all_chains = set()
    for info in subunit_info.values():
        all_chains.update(info["chain_names"])

    for chain_id in sorted(all_chains):
        if "error" in results.get(chain_id, {}):
            rows.append({
                "chain": chain_id,
                "is_splitted": None,
                "subunits": None,
                "total_missing": None,
                "error": results[chain_id]["error"],
            })
            continue

        chain_report = results[chain_id]

        # Collect subunits for this chain
        subunits = list(chain_report.keys())
        is_splitted = len(subunits) > 1

        # Count total missing residues across subunits
        total_missing = sum(len(info["missing"]) for info in chain_report.values())

        rows.append({
            "chain": chain_id,
            "is_splitted": is_splitted,
            "subunits": ",".join(subunits),
            "total_missing": total_missing,
        })

    df = pd.DataFrame(rows)
    return df
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font

def export_summary_excel_rich(results, subunit_info, filename="summary.xlsx"):
    rows = []
    all_chains = set()

    # Collect all chain names
    for info in subunit_info.values():
        all_chains.update(info["chain_names"])

    # Process each chain
    for chain_id in sorted(all_chains):
        if "error" in results.get(chain_id, {}):
            rows.append([chain_id, None, [], [], None, results[chain_id]["error"]])
            continue

        chain_report = results[chain_id]
        subunits_display = []
        missing_flags = []
        total_missing = 0

        for subunit, info in chain_report.items():
            total_missing += len(info["missing"])
            subunits_display.append(subunit)
            missing_flags.append(len(info["missing"]) > 0)

        is_splitted = len(chain_report) > 1
        rows.append([chain_id, is_splitted, subunits_display, missing_flags, total_missing, None])

    # Write with XlsxWriter
    with pd.ExcelWriter(filename, engine="xlsxwriter") as writer:
        df = pd.DataFrame(rows, columns=["chain", "is_splitted", "subunits", "missing_flags", "total_missing", "error"])

        # Create a display version with proper subunit formatting
        display_rows = []
        for _, row in df.iterrows():
            display_row = row.copy()
            # Convert single-item subunit lists to strings
            if isinstance(row["subunits"], list) and len(row["subunits"]) == 1:
                display_row["subunits"] = row["subunits"][0]
            display_rows.append(display_row)

        display_df = pd.DataFrame(display_rows).drop(columns=["missing_flags"])
        display_df.to_excel(writer, index=False, sheet_name="Summary")

        workbook = writer.book
        worksheet = writer.sheets["Summary"]

        # Create formats
        red_format = workbook.add_format({"color": "red"})
        black_format = workbook.add_format({"color": "black"})

        # Find "subunits" column index in the display DataFrame
        subunits_col = list(display_df.columns).index("subunits")

        # Apply rich formatting to subunits column
        for row_idx, (subs, flags) in enumerate(zip(df["subunits"], df["missing_flags"]), start=1):
            if not subs:  # skip rows with empty subunits (usually error rows)
                continue

            # Handle single subunit case vs multiple subunits
            if len(subs) == 1:
                # Single subunit - just apply formatting based on the single flag
                flag = flags[0] if flags else False
                if flag:
                    worksheet.write_rich_string(row_idx, subunits_col, red_format, subs[0])
                else:
                    worksheet.write_rich_string(row_idx, subunits_col, black_format, subs[0])
            else:
                # Multiple subunits - build rich string with separators
                rich_parts = []
                for i, (sub, flag) in enumerate(zip(subs, flags)):
                    if i > 0:
                        rich_parts.append(", ")  # add separator

                    # Apply appropriate formatting based on missing flag
                    if flag:
                        rich_parts.extend([red_format, sub])
                    else:
                        rich_parts.extend([black_format, sub])

                # Write rich string to the correct row (row_idx accounts for header)
                worksheet.write_rich_string(row_idx, subunits_col, *rich_parts)

    print(f"Excel summary with rich formatting saved to {filename}")


def get_pdb_long_jumps(pdb_file_path, chain_id, output_file=None):
    """
    Get a table of long jumps between consecutive Cα atoms from a PDB file.

    Parameters:
    -----------
    pdb_file_path : str
        Path to the PDB structure file
    chain_id : str, default 'J'
        Chain identifier to analyze
    distance_threshold : float, default 6.0
        Distance threshold in Angstroms to identify long jumps (typical CA-CA ~3.8 Å)
    output_file : str, optional
        Path to save the results table (CSV format). If None, no file is saved.

    Returns:
    --------
    pd.DataFrame: Table with columns ['Residue i', 'Residue i+1', 'Distance (Å)']
    """

    # Initialize PDB parser
    parser = PDBParser(QUIET=True)
    eps = 20.0

    try:
        # Parse the PDB file
        structure = parser.get_structure('protein', pdb_file_path)
        model = structure[0]
        chain = model[chain_id]

        # Extract Cα atoms and residues
        ca_atoms = []
        residues = []

        for residue in chain:
            # Skip hetero atoms and water molecules
            if residue.id[0] == ' ' and 'CA' in residue:
                ca_atoms.append(residue['CA'].get_coord())
                residues.append(residue)

        # Compute distances between consecutive Cα atoms
        long_jumps_data = []
        for i in range(1, len(ca_atoms)):
            residue_i = residues[i - 1].get_id()[1]
            residue_i_plus_1 = residues[i].get_id()[1]
            seq_dist = residue_i_plus_1 - residue_i
            distance = np.linalg.norm(ca_atoms[i] - ca_atoms[i - 1])
            distance_threshold = (seq_dist * 2) + eps

            # Check if this is a long jump
            if distance > distance_threshold:
                long_jumps_data.append({
                    'Residue i': residue_i,
                    'Residue i+1': residue_i_plus_1,
                    'Distance (Å)': round(distance, 2),
                    'Threshold(Å)': distance_threshold
                })

        # Create DataFrame
        long_jumps_df = pd.DataFrame(long_jumps_data)
        long_jumps_df = long_jumps_df.round(2)
        # Save to file if requested
        if output_file:
            long_jumps_df.to_csv(output_file, index=False)
            print(f"Long jumps table saved to: {output_file}")

        print(f"Found {len(long_jumps_df)} long jumps (>{distance_threshold} Å) in chain {chain_id}")

        return long_jumps_df

    except Exception as e:
        print(f"Error processing PDB file: {str(e)}")
        return pd.DataFrame()


if __name__ == "__main__":
    structure_file = "7e8t/cb_34_output_0.pdb"
    complex_name = "7e8t"
    #structure_file = "7e8t.pdb"
    subunit_json_file = "7e8t/subunits_info_after.json"
    #
    # with open(subunit_json_file, "r") as f:
    #     subunit_info = json.load(f)
    #
    # results = analyze_all_chains(structure_file, subunit_info)
    # # Build and display summary table
    # # summary_df = build_summary_table(results, subunit_info)
    # # print("\n=== Summary Table ===")
    # # print(summary_df.to_string(index=False))
    # # summary_df.to_csv(f"{complex_name}_summary.csv", index=False)
    # #
    # # # Print results
    # # for chain_id, chain_report in results.items():
    # #     print(f"\n=== Chain {chain_id} ===")
    # #     if "error" in chain_report:
    # #         print(" ", chain_report["error"])
    # #         continue
    # #     for subunit, info in chain_report.items():
    # #         if info["missing"]:
    # #             print(f" Subunit {subunit} (expected {info['expected_range']}):")
    # #             print(f"   - Missing {len(info['missing'])} residues "
    # #                   f"({info['missing'][0]}–{info['missing'][-1]})")
    # #         if info["mismatches"]:
    # #             print(f" Subunit {subunit} (expected {info['expected_range']}):")
    # #             for resnum, exp, obs in info["mismatches"]:
    # #                 print(f"   - Mismatch at {resnum}: expected {exp}, found {obs}")
    # #         # if not info["missing"] and not info["mismatches"]:
    # #         #     print("   - OK (all residues match)")
    # #     observed_chain = extract_chain_sequence(structure_file, chain_id)
    # #     print(f"Chain {chain_id} Structure len: {len(observed_chain)}, seq: {observed_chain}")
    # #     expected_seq, mapping = reconstruct_expected_chain_sequence(subunit_info, chain_id)
    # #     print(f"Chain {chain_id} Subunits Info len: {len(expected_seq)}, seq: {expected_seq}")
    # #     print(f"Done processing Chain {chain_id}")
    #
    # #export_summary_excel(results, subunit_info, filename=f"{complex_name}_summary.xlsx")
    # export_summary_excel_rich(results, subunit_info, filename=f"{complex_name}_summary.xlsx")
    chain_id = 'J'
    # Save to CSV file
    long_jumps_table = get_pdb_long_jumps(structure_file,chain_id, output_file=f'{complex_name}_long_jumps.csv')
    print(long_jumps_table)