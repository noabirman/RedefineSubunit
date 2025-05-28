import os
import json

def check_msa_against_subunits(base_dir):
    for complex_name in os.listdir(base_dir):
        complex_path = os.path.join(base_dir, complex_name)
        if not os.path.isdir(complex_path):
            continue

        subunits_info_path = os.path.join(complex_path, "subunits_info.json")
        msa_inputs_path = os.path.join(complex_path, "msa_inputs")
        mapping_path = os.path.join(complex_path, "chain_id_mapping.json")

        if not (os.path.isfile(subunits_info_path) and
                os.path.isdir(msa_inputs_path) and
                os.path.isfile(mapping_path)):
            print(f"Skipping {complex_name}: missing files")
            continue

        try:
            with open(subunits_info_path, 'r') as f:
                subunits_info = json.load(f)
            with open(mapping_path, 'r') as f:
                mapping = json.load(f)
        except Exception as e:
            print(f"Error reading files for {complex_name}: {e}")
            continue

        # Build reference chain sequences
        chain_ref_seqs = {}
        for sub in subunits_info.values():
            for chain in sub["chain_names"]:
                chain_ref_seqs.setdefault(chain, "")
                chain_ref_seqs[chain] += sub["sequence"]

        # Only initialize coverage for chain_ids actually used in the mapping
        used_chains = set(entry["chain_id"] for entry in mapping.values())
        chain_coverage = {chain: [False] * len(seq) for chain, seq in chain_ref_seqs.items() if chain in used_chains}

        for msa_file in os.listdir(msa_inputs_path):
            if not msa_file.endswith(".json"):
                continue

            msa_path = os.path.join(msa_inputs_path, msa_file)
            try:
                with open(msa_path, 'r') as f:
                    msa_data = json.load(f)
                msa_id = msa_data["name"]
                msa_seq = msa_data["sequences"][0]["protein"]["sequence"]
            except Exception as e:
                print(f"{complex_name}: Failed to parse {msa_file}: {e}")
                continue

            if msa_id not in mapping:
                print(f"{complex_name}: MSA ID {msa_id} not in mapping.")
                continue

            entry = mapping[msa_id]
            chain_id = entry["chain_id"]
            start, end = entry["start"], entry["end"]

            if chain_id not in chain_ref_seqs:
                print(f"{complex_name}: Chain {chain_id} from mapping not found in subunits_info.")
                continue

            ref_seq = chain_ref_seqs[chain_id]
            expected = ref_seq[start - 1:end]

            if msa_seq != expected:
                print(f"{complex_name}: MSA {msa_id} does not match reference {chain_id} {start}-{end}")
            else:
                for i in range(start - 1, end):
                    chain_coverage[chain_id][i] = True

        # Check for gaps
        for chain, covered in chain_coverage.items():
            if not all(covered):
                gaps = []
                start = None
                for i, flag in enumerate(covered):
                    if not flag and start is None:
                        start = i
                    elif flag and start is not None:
                        gaps.append((start + 1, i))
                        start = None
                if start is not None:
                    gaps.append((start + 1, len(covered)))
                print(f"{complex_name}: Chain {chain} missing coverage at {gaps}")


import itertools
from Bio.PDB.MMCIFParser import MMCIFParser

def get_chain_sequence(cif_path, chain_id):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", cif_path)
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                seq = ""
                for res in chain:
                    if res.id[0] == " " and "CA" in res:
                        resname = res.get_resname()
                        one_letter = residue_3_to_1.get(resname.upper(), 'X')
                        seq += one_letter
                return seq
    return None

# Standard 3-letter to 1-letter mapping
residue_3_to_1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def check_af_pair_completeness_and_integrity(base_dir):
    for complex_name in os.listdir(base_dir):
        complex_path = os.path.join(base_dir, complex_name)
        msa_path = os.path.join(complex_path, "msa_inputs")
        mapping_path = os.path.join(complex_path, "chain_id_mapping.json")
        af_pairs_path = os.path.join(complex_path, "af_pairs")
        subunits_info_path = os.path.join(complex_path, "subunits_info.json")  # ✅ NEW

        if not (os.path.isdir(msa_path) and os.path.isfile(mapping_path)
                and os.path.isdir(af_pairs_path) and os.path.isfile(subunits_info_path)):
            print(f"Skipping {complex_name}: missing required files or folders.")
            continue

        try:
            with open(mapping_path, 'r') as f:
                mapping = json.load(f)
            with open(subunits_info_path, 'r') as f:
                subunits_info = json.load(f)
        except Exception as e:
            print(f"{complex_name}: Failed to load mapping or subunits_info: {e}")
            continue

        # ✅ Detect which chain_ids came from multi-chain subunits (allow self-pairs only for them)
        self_pair_allowed_chain_ids = set()
        for subunit in subunits_info.values():
            if len(subunit.get("chain_names", [])) > 1:
                for msa_id, entry in mapping.items():
                    if set(subunit["chain_names"]) & {entry["chain_id"]}:
                        self_pair_allowed_chain_ids.add(entry["chain_id"].upper())

        # Load MSA sequences and their associated chain_ids
        msa_seqs_by_chain_id = {}
        for msa_file in os.listdir(msa_path):
            if not msa_file.endswith(".json"):
                continue
            with open(os.path.join(msa_path, msa_file)) as f:
                msa_data = json.load(f)
                msa_id = msa_data["name"]
                msa_seq = msa_data["sequences"][0]["protein"]["sequence"]
                if msa_id not in mapping:
                    print(f"{complex_name}: MSA ID {msa_id} missing from mapping.")
                    continue
                chain_id = mapping[msa_id]["chain_id"].upper()
                msa_seqs_by_chain_id[chain_id] = msa_seq

        chain_ids = sorted(msa_seqs_by_chain_id.keys())

        # ✅ Generate all unordered normalized expected pairs
        expected_pairs = set()
        for a, b in itertools.combinations(chain_ids, 2):
            pair = f"{min(a, b).lower()}_{max(a, b).lower()}"
            expected_pairs.add(pair)
        # Add valid self-pairs (only if allowed)
        for c in chain_ids:
            if c in self_pair_allowed_chain_ids:
                expected_pairs.add(f"{c.lower()}_{c.lower()}")

        found_pairs = set()
        mismatched_sequences = []
        extra_pairs = []

        # Validate af_pairs contents
        for pair_dir in os.listdir(af_pairs_path):
            pair_path = os.path.join(af_pairs_path, pair_dir)
            cif_path = os.path.join(pair_path, f"{pair_dir}_model.cif")

            if not os.path.isfile(cif_path):
                continue

            parts = pair_dir.split('_')
            if len(parts) != 2:
                print(f"{complex_name}: Invalid pair folder name {pair_dir}")
                continue

            a_id, b_id = parts[0].upper(), parts[1].upper()
            # Normalize actual found pair name
            normalized_pair = f"{min(a_id, b_id).lower()}_{max(a_id, b_id).lower()}"
            found_pairs.add(normalized_pair)

            if normalized_pair not in expected_pairs:
                extra_pairs.append(pair_dir)

            # Get CIF chain sequences
            seq_a = get_chain_sequence(cif_path, 'A')
            seq_b = get_chain_sequence(cif_path, 'B')
            expected_a = msa_seqs_by_chain_id.get(a_id)
            expected_b = msa_seqs_by_chain_id.get(b_id)

            mismatch = []
            if expected_a and expected_b:
                match_direct = (seq_a == expected_a and seq_b == expected_b)
                match_swapped = (seq_a == expected_b and seq_b == expected_a)
                if not (match_direct or match_swapped):
                    mismatch.append(f"Expected {a_id}/{b_id}, but sequences mismatched")
            else:
                if not expected_a:
                    mismatch.append(f"{a_id} missing in msa inputs")
                if not expected_b:
                    mismatch.append(f"{b_id} missing in msa inputs")

            if mismatch:
                mismatched_sequences.append((pair_dir, mismatch))

        # Report
        missing_pairs = sorted(expected_pairs - found_pairs)
        if missing_pairs:
            print(f"{complex_name}: Missing AF pairs: {missing_pairs}")
        if extra_pairs:
            print(f"{complex_name}: Extra/unexpected AF pairs: {extra_pairs}")
        if mismatched_sequences:
            for pair, mismatches in mismatched_sequences:
                print(f"{complex_name}: {pair}_model.cif sequence mismatch: {', '.join(mismatches)}")


if __name__ == "__main__":
    #check_msa_against_subunits("/cs/labs/dina/tsori/af3_example/complexes/splited_chains")
    check_af_pair_completeness_and_integrity("/cs/labs/dina/tsori/af3_example/complexes/DONE_MSA2")


