import pandas as pd
import os
from glob import glob
import json
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
import subprocess as sp

# --- CONFIGURATION ---
# Path to the big databases on the server
database_dir = "/cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/databases"
output_folder_name = "msa_lib"

# Define source database paths
pdb_seqres = f'{database_dir}/pdb_seqres_2022_09_28.fasta'
uniprot = f'{database_dir}/uniprot_all_2021_04.fa'
uniref = f'{database_dir}/uniref90_2022_05.fa'
mgy = f'{database_dir}/mgy_clusters_2022_05.fa'
bfd = f'{database_dir}/bfd-first_non_consensus_sequences.fasta'

MSA_refs = [bfd, mgy, uniprot, uniref, pdb_seqres]

# Global dictionary for headers
A_MSA_headers = dict()


def update_headers_one_seq(af3_json):
    """
    Parses an AlphaFold3 JSON to extract all protein IDs (templates and MSA hits).
    """
    try:
        data = json.load(open(af3_json, "r"))
        # Check if 'sequences' exists to avoid errors on empty/malformed JSONs
        if 'sequences' not in data:
            return

        protein_data = data['sequences'][0]['protein']

        # Parse Unpaired MSA
        A_uMSA = protein_data['unpairedMsa'].split("\n")[2:-1]
        for header_line in range(0, len(A_uMSA), 2):
            header = A_uMSA[header_line].split("/")[0][1:]
            seq = A_uMSA[header_line + 1].replace("-", "")
            A_MSA_headers[header] = seq

        # Parse Paired MSA
        A_pMSA = protein_data['pairedMsa'].split("\n")[2:-1]
        for header_line in range(0, len(A_pMSA), 2):
            header = A_pMSA[header_line].split("/")[0][1:]
            seq = A_pMSA[header_line + 1].replace("-", "")
            if header not in A_MSA_headers:
                A_MSA_headers[header] = seq

        # Parse Templates (PDB IDs)
        A_templates = protein_data.get('templates', [])
        for template in A_templates:
            line = template['mmcif'].split("\n")[0]
            pdb = (line.split("_")[1]).lower()
            if pdb not in A_MSA_headers:
                A_MSA_headers[pdb] = ""

    except Exception as e:
        print(f"Error reading {af3_json}: {e}")


def write_db_one(new_database_name, database):
    """
    Reads the huge database and writes valid entries to the msa_lib folder.
    """
    # Create the output path inside msa_lib
    output_path = os.path.join(output_folder_name, new_database_name)

    print(f"Writing {output_path}...")
    with open(output_path, 'w+') as fout:
        for record in SeqIO.parse(database, format='fasta'):
            # Normalize ID format
            if "pdb" in database:
                _id = str(record.id).split("_")[0]
            else:
                _id = str(record.id).split(" ")[0]

            # Write if in whitelist
            if _id in A_MSA_headers:
                SeqIO.write(record, fout, format='fasta')


if __name__ == '__main__':
    # 1. Create the msa_lib directory
    os.makedirs(output_folder_name, exist_ok=True)

    # 2. Collect Headers from all JSONs in current folder
    predictions = [file for file in os.listdir(".") if file.endswith("_data.json")]
    target_name = os.path.basename(os.getcwd())

    print(f"Found {len(predictions)} JSON files. Collecting headers...")
    for data_json in predictions:
        update_headers_one_seq(data_json)

    print(f"Total unique headers to find: {len(A_MSA_headers)}")

    # 3. Parallel processing to write filtered databases
    # Arguments: (Desired Filename, Source Database Path)
    tasks = [(f"{target_name}_{os.path.basename(seq_db)}", seq_db) for seq_db in MSA_refs]

    with Pool(min(cpu_count(), len(MSA_refs))) as p:
        p.starmap(write_db_one, tasks)

    print("Finished filtering proteins.")

    # 4. Create dummy RNA files inside msa_lib
    rfam = "rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta"
    rnacentral = "rnacentral_active_seq_id_90_cov_80_linclust.fasta"
    nt_rna = "nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta"
    empty_msa_refs = [rfam, rnacentral, nt_rna]

    print("Creating dummy RNA files...")
    for empty_msa_ref in empty_msa_refs:
        # Create empty file inside msa_lib
        sp.run(f"touch {output_folder_name}/{empty_msa_ref}", shell=True)

    # 5. Create symbolic link for mmcif inside msa_lib
    print("Linking mmcif directory...")
    # Check if link exists to avoid error
    if not os.path.exists(f"{output_folder_name}/mmcif_files"):
        sp.run(f"ln -s {database_dir}/mmcif_files {output_folder_name}/mmcif_files", shell=True)

    print("Done. All files are in:", os.path.abspath(output_folder_name))