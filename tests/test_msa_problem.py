import json
import sys

def compare_sequences(msa_file, data_file, chain_id):
    # Load AF3 input JSON
    with open(data_file, "r") as f:
        data = json.load(f)

    query_seq = None
    for seq in data["sequences"]:
        if seq["protein"]["id"] == chain_id:
            query_seq = seq["protein"]["sequence"]
            break

    if query_seq is None:
        print(f"Chain {chain_id} not found in {_data.json}")
        return

    # Load MSA JSON
    with open(msa_file, "r") as f:
        msa_data = json.load(f)

    if "msa" not in msa_data:
        print(f"No 'msa' key found in {msa_file}")
        return

    msa = msa_data["msa"]

    # Determine if paired or unpaired
    if "paired_msa" in msa:
        print("MSA type: paired")
        # assume chain_id is second in the pair -> row_b
        first_seq = msa["paired_msa"][0]["row_b"]["sequence"]
    elif "unpaired_msa" in msa:
        print("MSA type: unpaired")
        first_seq = msa["unpaired_msa"][0]["sequence"]
    else:
        print("No paired or unpaired MSA found")
        return

    # Compare sequences
    if first_seq == query_seq:
        print("✅ Sequences match exactly!")
    else:
        print("❌ Mismatch detected!")
        # find position of first mismatch
        for i, (q, m) in enumerate(zip(query_seq, first_seq)):
            if q != m:
                print(f"First mismatch at position {i}: query='{q}' msa='{m}'")
                break
        print(f"Query seq (start): {query_seq[:50]}...")
        print(f"MSA   seq (start): {first_seq[:50]}...")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python test_msa_match.py <msa_json> <data_json> <chain_id>")
        sys.exit(1)
    msa_file = sys.argv[1]
    data_file = sys.argv[2]
    chain_id = sys.argv[3]

    compare_sequences(msa_file, data_file, chain_id)
