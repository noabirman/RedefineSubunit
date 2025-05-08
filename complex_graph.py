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


# def extract_sequence_with_seqio(structure_path,af_version: int):
#     """
#     Extracts the sequence from an mmCIF/PDB file using Bio.SeqIO.
#
#     Args:
#         structure_path (str): Path to the mmCIF/PDB file.
#         af_version if 2 then PDB and if 3 cif
#
#     Returns:
#         str: The amino acid sequence as a single-letter code string.
#     """
#     format = {'2':"pdb-atom", '3':"cif-atom"}
#     sequences = []
#     for record in SeqIO.parse(structure_path, format[af_version]):
#         sequences.append(str(record.seq))
#         print(record.id)
#     return ''.join(sequences)

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


def find_high_confidence_regions(plddt_array, confidence_threshold=40, gap_threshold=3):
    """
    Finds ranges of high-confidence regions in a plDDT array while preserving the order.

    Args:
        plddt_array (np.ndarray): Array containing plDDT values for each residue.
        confidence_threshold (float): Minimum confidence value to include a residue.
        gap_threshold (int): Maximum allowed gap between indices in the same region.

    Returns:
        list[tuple[int, int]]: A list of tuples, each representing the start and end indices of a high-confidence region.
    """
    indices = np.where(plddt_array > confidence_threshold)[0]

    if len(indices) == 0:
        return []

    regions = []
    start_index = indices[0]

    for i in range(1, len(indices)):
        if indices[i] - indices[i - 1] > gap_threshold:
            if indices[i - 1] - start_index >= 5:
                regions.append((int(start_index), int(indices[i - 1])))
            start_index = indices[i]

    # Handle the last region (fix: use indices[-1])
    if indices[-1] - start_index >= 5:
        regions.append((int(start_index), int(indices[-1])))

    return regions



def extract_subunit_info(indexs: List[Tuple[int, int]], token_chain_ids: List[str], full_seq: str) -> List[SubunitInfo]:
    """
    Build a list of SubunitInfo ensuring different chains are treated as separate subunits.

    Args:
        indexs (List[Tuple[int, int]]): List containing the index of each subunit.
        token_chain_ids (List[str]): Indicate the chain IDs of each subunit.
        full_seq (str): The residue sequence of the current protein.

    Returns:
        list[SubunitInfo]: A list of SubunitInfo objects.
    """
    subunit_infos = []
    chain_occ_counter = Counter()  # Track occurrences for each chain separately

    for start, end in indexs:
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
                chain_names=[chain_id],  # Only one chain per subunit
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
        threshold (int): adding edge iff pae of the link between two vertices > threshold.

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
    """Renames chain identifiers based on filename components.

    Args:
        data_path (Path): Path to the data file (used for extracting chain names).
        token_chain_ids (list[str]): Original chain identifiers.

    Returns:
        list[str]: Updated chain identifiers.
    """

    # Extract filename without extension
    #filename = data_path.stem  # "cdf_afg_ajs_confidences" → "cdf_afg_ajs"
    filename = os.path.basename(data_path)  # Extract 'cdf_ddf' #todo: assume here there is only 2 chains
    # Split by '_' to get chain names
    parts = filename.split('_')

    # Ensure we have enough parts to map chains
    if len(parts) < 2:
        raise ValueError(f"Unexpected filename format: {filename}. Expected at least 2 parts.")

    # Generate a mapping dynamically (A → first part, B → second, etc.)
    # In case of chain with itself (e_e ) change to E1, E2
    if parts [0] == parts[1]: # chain with itself case
        replacement_dict = {chr(65 + i): part.upper() + "_"+ str (i+1) for i, part in enumerate(parts)}
    else: #regular case
        replacement_dict = {chr(65 + i): part.upper() for i, part in enumerate(parts)}

    # Update token_chain_ids based on mapping
    return [replacement_dict.get(chain, chain) for chain in token_chain_ids]


def graph(structure_path: str, data_path:str, af_version: str)->nx.Graph: # het
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
    groups_indexs = find_high_confidence_regions(plddt_array, confidence_threshold=40)

    token_chain_ids_updated = rename_chains_from_file(data_path, token_chain_ids)

    subunits_info = extract_subunit_info(groups_indexs, token_chain_ids_updated, full_seq)
    G = nx.Graph()
    edges = find_edges(subunits_info, pae_as_arr, threshold=15)
    # index per chain instead of per full seq (by doing token_res_ids[subunit.start])
    updated_subunits = [
        dataclasses.replace(subunit, start=token_res_ids[subunit.start], end=token_res_ids[subunit.end])
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
