import json


def json_to_fasta(json_file:str): #should be path/synapse.json
    # Load your JSON file
    with open(json_file) as f:
        data = json.load(f)

    # Write to a FASTA file
    with open("sequences.fasta", "w") as fasta:
        for entry in data["sequences"]:
            pid = entry["protein"]["id"]
            seq = entry["protein"]["sequence"]
            fasta.write(f">{pid}\n{seq}\n")


if __name__ == '__main__':
    json_to_fasta("/cs/labs/dina/tsori/af3_example/msa_inputs/synapse.json")

