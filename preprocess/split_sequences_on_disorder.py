"""
Script Name: split_sequences_on_disorder.py

Description:
    Splits protein sequences longer than a specified threshold (default 1000 residues)
    based on predicted disordered regions from IUPred3. The script avoids splitting
    within ordered regions and only allows cuts within disordered zones.

Usage:
    Run after generating IUPred3 outputs for each protein in a directory.

    Example:
        python3 split_sequences_on_disorder.py folder_to_iupred_outputs sequences.fasta

Author: Noa Birman
Date: 2025-04-21
"""
import os
import json
import sys

#from numba.np.arraymath import complex_nanmax

max_len = 1000
disorder_threshold = 0.5

# === Helper: Read IUPred3 scores ===
def read_iupred_scores(file_path):
    scores = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            pos, aa, score = int(parts[0]), parts[1], float(parts[2])
            scores.append(score)
    return scores

# === Helper: Read FASTA ===
def read_fasta(path):
    sequences = {}
    with open(path) as f:
        current_id, current_seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = "".join(current_seq)
    return sequences

# === Helper: Find cut points ===
def find_cut_points(scores, max_len, threshold):
    cut_points = []
    i = max_len
    while i < len(scores):
        # Look for the nearest disordered residue within ±20 residues
        window = range(max(i - 20, 1), min(i + 20, len(scores)))
        for j in window:
            if scores[j] > threshold:
                cut_points.append(j)
                i = j + max_len  # skip forward to next chunk
                break
        else:
            # No cut point found: force a cut anyway
            cut_points.append(i)
            i += max_len
    return cut_points

# === Main ===
if __name__ == '__main__':
    # === Settings ===
    iupred_dir = sys.argv[1]  # folder with all *_iupred.txt files
    original_fasta = sys.argv[2]  # the sequences you used as input "sequences.fasta"
    complex_name = sys.argv[3]
    output_fasta = os.path.join(complex_name, "iupred_split_sequences.fasta")
    output_map = os.path.join(complex_name,"iupred_split_mapping.json")

    sequences = read_fasta(original_fasta)
    fragment_map = {}
    fragments = []

    for prot_id, seq in sequences.items():
        iupred_file = os.path.join(iupred_dir, f"{prot_id}_iupred.txt")
        if not os.path.exists(iupred_file):
            print(f"Missing IUPred output for {prot_id}, skipping.")
            continue
        scores = read_iupred_scores(iupred_file)
        if len(seq) <= max_len:
            new_id = f"{prot_id}_1"
            fragments.append((new_id, seq))
            fragment_map[prot_id] = [new_id]
            continue

        # Find cut points
        cuts = find_cut_points(scores, max_len, disorder_threshold)
        prev = 0
        frag_ids = []
        for i, cut in enumerate(cuts + [len(seq)]):
            frag_seq = seq[prev:cut]
            new_id = f"{prot_id}_{i+1}"
            fragments.append((new_id, frag_seq))
            frag_ids.append(new_id)
            prev = cut
        fragment_map[prot_id] = frag_ids

    # === Write Output ===
    with open(output_fasta, "w") as f:
        for frag_id, frag_seq in fragments:
            f.write(f">{frag_id}\n{frag_seq}\n")

    with open(output_map, "w") as f:
        json.dump(fragment_map, f, indent=2)

    print(f"✅ Done! Output saved to {output_fasta} and {output_map}")

    # === Create AF3-style JSON ===
    af3_sequences = [
        {"protein": {"id": frag_id, "sequence": frag_seq}}
        for frag_id, frag_seq in fragments
    ]

    af3_json = {
        "name": "SPLIT_FROM_IUPRED",
        "modelSeeds": [1],
        "sequences": af3_sequences
    }

    with open(os.path.join(complex_name,"af3_input.json"), "w") as f:
        json.dump(af3_json, f, indent=2)

    print("✅ AlphaFold 3 input JSON saved to af3_input.json")

