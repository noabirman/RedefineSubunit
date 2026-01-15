import numpy as np
import json
import sys
import os
import itertools
from collections import defaultdict, Counter
from typing import List, Tuple
from Bio.PDB import MMCIFParser
import Bio.SeqUtils
import Bio.PDB, Bio.PDB.Residue
from Bio import SeqIO, BiopythonParserWarning
import dataclasses
import networkx as nx
from pathlib import Path
import warnings

SubunitName = str
@dataclasses.dataclass
class SubunitInfo:
    name: SubunitName
    chain_names: List[str]
    start: int  # was before indexs: Tuple[int, int]
    end: int
    sequence: str

def parse_ipsae_scores(ipsae_path: str) -> np.ndarray:
    """
    Parses the ipSAE _byres.txt file.
    Assumes the last column is the ipSAE score (psae_d0res).
    """
    scores = []
    with open(ipsae_path, 'r') as f:
        next(f)  # Skip the header line
        for line in f:
            if not line.strip(): continue
            parts = line.split()
            # The last column in the _byres.txt file is the relevant ipSAE score
            try:
                scores.append(float(parts[-1]))
            except ValueError:
                continue
    return np.array(scores)

def extract_sequence_with_seqio(structure_path, af_version: int):
    """
    Extracts the sequence from an mmCIF/PDB file using Bio.SeqIO.

    Args:
        structure_path (str): Path to the mmCIF/PDB file.
        af_version (int): If 2, parse as PDB; if 3, parse as CIF (AlphaFold v3).

    Returns:
        str: The amino acid sequence as a single-letter code string.
    """
    format_map = {'2': "pdb-atom", '3': "cif-atom"}
    format = format_map.get(af_version)

    sequences = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonParserWarning)
        for record in SeqIO.parse(structure_path, format):
            sequences.append(str(record.seq))
            #print(record.id or "No ID")

    return ''.join(sequences)


def find_high_confidence_regions(plddt_array, chain_ids, confidence_threshold=40, gap_threshold=3):
    """
    Finds ranges of high-confidence regions in a plDDT array while preserving the order
    and ensuring regions are within the same chain.

    Args:
        plddt_array (np.ndarray): Array containing plDDT values for each residue.
        chain_ids (list[str]): List of chain IDs corresponding to each residue.
        confidence_threshold (float): Minimum confidence value to include a residue.
        gap_threshold (int): Maximum allowed gap between indices in the same region.

    Returns:
        list[tuple[int, int]]: A list of tuples, each representing the start and end 
                               indices of a high-confidence region.
    """
    regions = []
    unique_chains = set(chain_ids)

    for chain in unique_chains:
        indices = [i for i, (score, c) in enumerate(zip(plddt_array, chain_ids)) if score > confidence_threshold and c == chain]

        if not indices:
            continue

        start_index = indices[0]

        for i in range(1, len(indices)):
            if indices[i] - indices[i - 1] > gap_threshold:
                if indices[i - 1] - start_index >= 5:
                    regions.append((start_index, indices[i - 1]))
                start_index = indices[i]

        # Handle the last region
        if indices[-1] - start_index >= 5:
            regions.append((start_index, indices[-1]))

    return regions

def extract_subunit_info(indices: List[Tuple[int, int]], token_chain_ids: List[str], full_seq: str) -> List[SubunitInfo]:
    """
    Build a list of SubunitInfo ensuring different chains are treated as separate subunits.

    Args:
        indices (List[Tuple[int, int]]): List containing the index of each subunit.
        token_chain_ids (List[str]): Indicate the chain IDs of each subunit.
        full_seq (str): The residue sequence of the current protein.

    Returns:
        list[SubunitInfo]: A list of SubunitInfo objects.
    """
    subunit_infos = []
    chain_occ_counter = Counter()  # Track occurrences for each chain separately

    for start, end in indices:
        # Find unique chain IDs in the segment (preserving order)
        chains_ids_in_node = list(dict.fromkeys(token_chain_ids[start:end + 1]))

        # Process each chain separately
        for chain_id in chains_ids_in_node:
            # Identify residues belonging to this chain in the given range
            chain_positions = [i for i in range(start, end + 1) if token_chain_ids[i] == chain_id]

            if not chain_positions:  # Skip if no positions found (shouldn't happen)
                continue

            chain_start, chain_end = chain_positions[0], chain_positions[-1]
            subunit_name = f"{chain_id}_{chain_occ_counter[chain_id] + 1}"  # Unique numbering per chain

            # Increment the occurrence counter for this chain
            chain_occ_counter[chain_id] += 1

            subunit_infos.append(SubunitInfo(
                name=subunit_name,
                chain_names=[chain_id.split('_')[0]],  # Only one chain per subunit
                start=chain_start,
                end=chain_end,
                sequence=full_seq[chain_start:chain_end + 1]
            ))

    return subunit_infos


def atom_plddt_to_res_plddt(structure, atom_plddts:List[float]):
    """
        Convert plddt list in AF3 case to be per residue instead of per atom, calculate the average plddt for each res.

        Args:
            structure: cif structure from AF
            atom_plddts (List[float]): plddt per atom

        Returns:
            List[float]: plddt per residue.
        """
    # Map atoms to residues
    residue_plddt_sum = defaultdict(float)
    residue_atom_count = defaultdict(int)

    atom_index = 0  # Track atom index for atom_plddts

    for model in structure:
        for chain in model:
            for residue in chain:
                residue_key = (chain.id, residue.id[1])  # Use (chain ID, residue number) as key
                for atom in residue:
                    if atom_index < len(atom_plddts):
                        residue_plddt_sum[residue_key] += atom_plddts[atom_index]
                        residue_atom_count[residue_key] += 1
                        atom_index += 1
                    else:
                        print(f"Warning: atom_index {atom_index} exceeds atom_plddts length.")
                        break

    # Calculate average plDDT for each residue
    average_residue_plddt = {
        key: residue_plddt_sum[key] / residue_atom_count[key]
        for key in residue_plddt_sum
    }
    # Assuming 'average_residue_plddt' and 'pae_as_arr' are already defined
    # Convert average_residue_plddt to a list in the correct order
    plddt_values = [average_residue_plddt[key] for key in sorted(average_residue_plddt)]
    return np.array(plddt_values)


def find_edges(subunits_info: List[SubunitInfo], pae_matrix: np.array, threshold: int = 15) -> list[tuple[str, str, float]]:
    """
    Find edges by pae_matrix and threshold. adding edge iff pae of the link between two vertices > threshold.

    Args:
        subunits_info (List[SubunitInfo]): vertices information.
        pae_matrix (np.array): PAE matrix.
        threshold (int): adding edge iff pae of the link between two vertices < threshold.

    Returns:
        List[tuple[str,str,float]]: edges list, each edge is a tuple (v1, v2, weight).
    """
    edges = []
    for subunit1, subunit2 in itertools.combinations(subunits_info, 2):
        pae_rect = pae_matrix[subunit1.start:subunit1.end, subunit2.start:subunit2.end]
        if pae_rect.size == 0: #todo:not sure if best practice (M is start=132 end =132 so the rect comes out empty)
            continue
        pae_score = np.mean(pae_rect)
        if pae_score < threshold:
            edges.append((subunit1.name, subunit2.name, float(pae_score)))

    return edges


def get_chain_ids_per_residue(structure):
    """
    Made for getting the token_chain_ids which require for extract_subunit_info() in AF2 case.
    token_chain_ids len equals to the sequance len (number of residues).

    Args:
        structure (int): structure from PDB file.

    Returns:
        List[char]: token_chain_ids will look like ['A','A',..,'A','B','B',..,'B',..,'E','E',..,'E']
    """
    chain_ids = []  # List to store chain IDs per residue

    for model in structure:  # Iterate through models (usually only 1)
        for chain in model:  # Iterate through chains
            for residue in chain:  # Iterate through residues
                if residue.id[0] == " ":  # Exclude heteroatoms (like water, ligands)
                    chain_ids.append(chain.id)  # Store chain ID per residue

    return chain_ids


def rename_chains_from_file(data_path: str, token_chain_ids: list[str]) -> list[str]:
    filename = os.path.basename(data_path)
    # Assumes filename structure is "Name1_Name2_..."
    target_names = filename.split('_')[:2]  # Adjust slice if you have more chains

    # 1. Detect chains in the order they appear in the JSON
    # This captures ['A', 'B'] OR ['B', 'I'] without assuming specific letters
    unique_json_chains = list(dict.fromkeys(token_chain_ids).keys())

    # 2. Safety Check: Do residue counts match filename parts?
    if len(unique_json_chains) != len(target_names):
        print(f"Warning: JSON has {len(unique_json_chains)} chains but filename implies {len(target_names)}.")
        # FALLBACK: If counts mismatch, trust the JSON and do not rename
        return token_chain_ids

    # 3. Create the map based on ORDER
    # 1st JSON Chain -> 1st Filename Part
    # 2nd JSON Chain -> 2nd Filename Part
    replacement_dict = {
        old: new.upper()
        for old, new in zip(unique_json_chains, target_names)
    }

    # 4. Apply
    return [replacement_dict.get(c, c) for c in token_chain_ids]


def token_to_residue_pae(pae_as_arr: np.ndarray, token_res_ids: list[int], token_chain_ids: list[str]) -> np.ndarray:
    """
    Optimized version: Reduces matrix in two linear passes (2N) instead of quadratic (N^2).
    """
    # 1. Map (chain_id, res_id) to token indices
    res_key_to_tokens = defaultdict(list)
    for token_idx, (res_id, chain_id) in enumerate(zip(token_res_ids, token_chain_ids)):
        res_key_to_tokens[(chain_id, res_id)].append(token_idx)

    # 2. Prepare sorted keys (residues)
    sorted_unique_keys = list(res_key_to_tokens.keys())
    n_res = len(sorted_unique_keys)
    n_tokens = pae_as_arr.shape[0]

    # 3. Pass 1: Collapse Rows (Average over the tokens belonging to residue i)
    # Output shape: (n_res, n_tokens)
    semi_reduced = np.zeros((n_res, n_tokens), dtype=np.float32)

    for i, key in enumerate(sorted_unique_keys):
        tokens_indices = res_key_to_tokens[key]
        # Take the mean of the rows for this residue
        semi_reduced[i, :] = np.mean(pae_as_arr[tokens_indices, :], axis=0)

    # 4. Pass 2: Collapse Columns (Average over the tokens belonging to residue j)
    # Output shape: (n_res, n_res)
    final_res_pae = np.zeros((n_res, n_res), dtype=np.float32)

    for j, key in enumerate(sorted_unique_keys):
        tokens_indices = res_key_to_tokens[key]
        # Take the mean of the columns for this residue (from the semi_reduced matrix)
        final_res_pae[:, j] = np.mean(semi_reduced[:, tokens_indices], axis=1)

    return final_res_pae


def graph(structure_path: str, data_path: str, af_version: str, ipsae_path: str = None) -> nx.Graph:
    # args: "fold_mll4_1100_end_rbbp5_wdr5_p53x2/fold_mll4_1100_end_rbbp5_wdr5_p53x2_model_0.cif" "fold_mll4_1100_end_rbbp5_wdr5_p53x2/fold_mll4_1100_end_rbbp5_wdr5_p53x2_full_data_0.json" 3
    # args: "example/cdf_ddf/cdf_ddf_model.cif" "example/cdf_ddf/cdf_ddf_confidences.json" 3
    with open(data_path, "r") as file:
        json_full_data = json.load(file)
    pae_as_arr = np.array(json_full_data['pae'])
    if af_version == '3':
        atom_plddts = json_full_data['atom_plddts']
        atom_chain_ids = json_full_data['atom_chain_ids']
        token_res_ids = json_full_data['token_res_ids']
        token_chain_ids = json_full_data['token_chain_ids']  # per res
        # Parse the CIF file
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", structure_path)
        plddt_array = atom_plddt_to_res_plddt(structure, atom_plddts)
    elif af_version == '2':
        # json data include ['max_pae', 'pae', 'plddt', 'ptm', 'iptm']
        # plddt per res and not per atom
        plddt_array = np.array(json_full_data['plddt'])  # todo: not sure if it keeps order this way
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("original_pdb", structure_path)
        token_chain_ids = get_chain_ids_per_residue(structure)
    full_seq = extract_sequence_with_seqio(structure_path,
                                           af_version)
    #print(full_seq)

    token_chain_ids_updated = rename_chains_from_file(data_path, token_chain_ids)

    # phosporilation bug fix:
    residue_chain_ids = []
    residue_numbers = []
    seen_residues = set()  # Stores tuples of (chain_id, res_id)

    for res_id, chain_id in zip(token_res_ids, token_chain_ids_updated):
        unique_key = (chain_id, res_id)
        if unique_key not in seen_residues:
            residue_chain_ids.append(chain_id)
            residue_numbers.append(res_id)
            seen_residues.add(unique_key)

    # === MODIFICATION START ===
    if ipsae_path and os.path.exists(ipsae_path):
        print(f"Cutting subunits based on ipSAE scores from: {os.path.basename(ipsae_path)}")
        # Use ipSAE scores.
        # Threshold > 0 means the residue contributes to the interaction.
        score_array = parse_ipsae_scores(ipsae_path)

        # Safety check: Ensure array lengths match
        if len(score_array) != len(residue_chain_ids):
            raise ValueError(
                f"Mismatch: ipSAE file has {len(score_array)} scores, but structure has {len(residue_chain_ids)} residues.")

        # Use threshold 0.0 (anything > 0 is part of the interface)
        groups_indices = find_high_confidence_regions(score_array, residue_chain_ids, confidence_threshold=0.0)
    else:
        print("Cutting subunits based on pLDDT scores (Default)")
        # Use existing pLDDT array with standard threshold (40)
        groups_indices = find_high_confidence_regions(plddt_array, residue_chain_ids, confidence_threshold=40)
    # === MODIFICATION END ===

    # Suppose pae_as_arr is n_tokens x n_tokens, token_res_ids contains repeated tokens
    residue_pae = token_to_residue_pae(pae_as_arr, token_res_ids, token_chain_ids_updated)
    # Now residue_pae.shape == (num_residues, num_residues

    subunits_info = extract_subunit_info(groups_indices, residue_chain_ids, full_seq)

    G = nx.Graph()
    edges = find_edges(subunits_info, residue_pae, threshold=15)
    # index per chain instead of per full seq (by doing token_res_ids[subunit.start])
    res_num = np.unique(token_res_ids) #addition for phosp
    updated_subunits = [
        dataclasses.replace(subunit, start=residue_numbers[subunit.start], end=residue_numbers[subunit.end])
        for subunit in subunits_info
    ]

    for subunit in updated_subunits:
            G.add_node(subunit.name, data=subunit)
    for e in edges: # e is (v1, v2, weight)
        G.add_edge(e[0], e[1], weight=e[2])
    return G


if __name__ == '__main__':
    if len(sys.argv) == 4:
        structure_path, data_path, af_version = os.path.abspath(sys.argv[1]),os.path.abspath(sys.argv[2]),sys.argv[3]
    else:
        print("usage: <script> structure_path data_path af_version")
    g = graph(structure_path, data_path, af_version)
    print (g)
    #plot_pae_plddt(pae_as_arr, plddt_array, nodes_as_req, edges, 'skip4_pae15_')
    #plot_pae_plddt2(pae_as_arr, plddt_array, nodes_as_req, edges, 'with_weights')
