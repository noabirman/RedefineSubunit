import os
from functools import partial
from glob import glob
import json
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from multiprocessing import Pool
import sys

from tqdm import tqdm

# Path to the databases
database_dir = "/cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/databases"

pdb_seqres = f'{database_dir}/pdb_seqres_2022_09_28.fasta'
uniprot = f'{database_dir}/uniprot_all_2021_04.fa'
uniref = f'{database_dir}/uniref90_2022_05.fa'
mgy = f'{database_dir}/mgy_clusters_2022_05.fa'
bfd = f'{database_dir}/bfd-first_non_consensus_sequences.fasta'

MSA_DBS = [bfd, mgy, uniprot, uniref, pdb_seqres]

EMPTY_MSA_DBS = ["rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta",
                 "rnacentral_active_seq_id_90_cov_80_linclust.fasta",
                 "nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta"
    ]

MMCIF_FOLDER = f'{database_dir}/mmcif_files'


def get_records_from_json(af3_json):
    data = json.load(open(af3_json, "r"))
    record_to_keep = set()

    for seq in data['sequences']:
        if "protein" not in seq:
            continue
        seq_data = seq['protein']
        unpaired_msa = seq_data['unpairedMsa'].split("\n")[2:-1]
        paired_msa = seq_data['pairedMsa'].split("\n")[2:-1]
        templates = seq_data['templates']

        for header_line_ind in range(0, len(unpaired_msa), 2):
            header = unpaired_msa[header_line_ind].split("/")[0][1:]
            # seq = unpaired_msa[header_line_ind + 1].replace("-", "")
            # A_MSA_headers[header] = seq
            record_to_keep.add(header)

        for header_line_ind in range(0, len(paired_msa), 2):
            header = paired_msa[header_line_ind].split("/")[0][1:]
            record_to_keep.add(header)

            for template in templates:
                line = template['mmcif'].split("\n")[0]
                pdb = (line.split("_")[1]).lower()
                record_to_keep.add(pdb)
    return record_to_keep


def process_db(database_path, records_to_keep, output_folder):
    database_basename = os.path.basename(database_path)
    print(f"Processing database: {database_basename}")
    new_database_name = os.path.join(output_folder, database_basename)
    total_records = 0
    with open(new_database_name, 'w+') as fout:
        for record in SeqIO.parse(database_path, format='fasta'):
            if "pdb" in database_path:
                _id = str(record.id).split("_")[0]
            else:
                _id = str(record.id).split(" ")[0]
            if _id in records_to_keep:
                SeqIO.write(record, fout, format='fasta')
                total_records += 1
    print(f"Wrote {new_database_name} with {total_records} records")


def main(jsons_folder, output_folder):
    jsons_folder = os.path.abspath(jsons_folder)
    output_folder = os.path.abspath(output_folder)

    # use glob to find all _data.json in all inner folders
    predictions = glob(os.path.join(jsons_folder, '**', '*_data.json'), recursive=True)
    os.makedirs(output_folder, exist_ok=True)

    # predictions = predictions[:20]

    print("Collecting records to keep...")
    with Pool() as pool:
        all_records_set = list(tqdm(pool.imap(get_records_from_json, predictions), total=len(predictions)))
    # all_records_set = []
    # for data_json in tqdm(predictions):
    #     records = get_records_from_json(data_json)
    #     all_records_set.append(records)

    cache_path = os.path.join(output_folder, 'records_to_keep_cache.json')
    json.dump([list(s) for s in all_records_set], open(cache_path, 'w'), indent=2)

    records_to_keep = set().union(*all_records_set)
    print(f"Total records to keep: {len(records_to_keep)}")

    record_counter = {}
    for record_set in all_records_set:
        for record in record_set:
            if record not in record_counter:
                record_counter[record] = 0
            record_counter[record] += 1
    for percentage in (0.1, 0.5, 0.9, 0.98):
        threshold = int(len(all_records_set) * percentage)
        num_records = sum(1 for count in record_counter.values() if count >= threshold)
        print(f"Records present in at least {percentage*100}% of jsons ({threshold} files): {num_records} ({num_records/len(record_counter)*100:.2f}%)")

    print("Writing new databases...")
    process_db_partial = partial(process_db, records_to_keep=records_to_keep, output_folder=output_folder)
    with Pool() as pool:
        list(tqdm(pool.imap(process_db_partial, MSA_DBS), total=len(MSA_DBS)))
    # for database_path in tqdm(MSA_DBS):
    #     process_db(database_path, records_to_keep, output_folder)

    print("Preparing the dataset needed files")
    for empty_db in EMPTY_MSA_DBS:
        new_empty_db_path = os.path.join(output_folder, empty_db)
        open(new_empty_db_path, 'w+').close()
        print(f"Created empty database file: {new_empty_db_path}")
    # link to mmcif folder
    new_mmcif_folder = os.path.join(output_folder, 'mmcif_files')
    if not os.path.exists(new_mmcif_folder):
        os.symlink(MMCIF_FOLDER, new_mmcif_folder)
        print(f"Linked mmcif_files folder to: {new_mmcif_folder}")

    print("Done.")


if __name__ == '__main__':
    assert len(sys.argv) == 3, "Usage: <script> jsons_folder output_folder"

    main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]))