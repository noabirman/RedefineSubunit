import os
import json
from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder

def extract_chain_sequences(pdb_file):
    """Return dict of {chain_id: sequence} from a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    ppb = PPBuilder()
    chain_seqs = {}
    for model in structure:
        for chain in model:
            seq = ""
            for pp in ppb.build_peptides(chain):
                seq += str(pp.get_sequence())
            if seq:
                chain_seqs[chain.id] = seq
    return chain_seqs

def search_sequence_in_pdbs(folder, query_seq):
    """Search all PDBs in folder for chains containing query_seq."""
    query_seq = query_seq.upper()
    matches = []
    for file in os.listdir(folder):
        if not file.endswith(".pdb"):
            continue
        pdb_path = os.path.join(folder, file)
        chain_seqs = extract_chain_sequences(pdb_path)
        for chain_id, seq in chain_seqs.items():
            if query_seq in seq:
                matches.append({
                    "pdb_file": file,
                    "chain_id": chain_id,
                    "match_position": seq.find(query_seq) + 1,
                    "chain_length": len(seq)
                })
    return matches

def search_subunits_in_pdb_folder(subunits_json, pdb_folder):
    with open(subunits_json) as f:
        subunits = json.load(f)

    results = {}
    for sub_name, info in subunits.items():
        seq = info["sequence"]
        print(f"🔎 Searching for subunit {sub_name} ({len(seq)} aa)...")
        matches = search_sequence_in_pdbs(pdb_folder, seq)
        results[sub_name] = matches
        if matches:
            for m in matches:
                print(f"  ✅ Found in {m['pdb_file']} chain {m['chain_id']} (pos {m['match_position']})")
        else:
            print("  ❌ No match found.")
        print()
    return results

if __name__ == "__main__":
    # --- edit these paths ---
    subunits_json = "/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2/8a5o_again/combfold/subunits_info.json"
    pdb_folder = "/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2/8a5o_again/combfold/models"
    output_file = "8a5o_old_missing_subunit_in_models_results.json"

    results = search_subunits_in_pdb_folder(subunits_json, pdb_folder)

    # Save all matches to a JSON
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n✅ Results saved to {output_file}")