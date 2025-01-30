import numpy as np
import json
from Bio.PDB import MMCIFParser
from collections import defaultdict
from typing import List, Tuple
from collections import Counter
import sys
import os
import Bio.SeqUtils
import Bio.PDB, Bio.PDB.Residue
from Bio import SeqIO
from Bio.PDB import MMCIFParser, PPBuilder
import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
import dataclasses

SubunitName = str


@dataclasses.dataclass
class SubunitInfo:
    name: SubunitName
    chain_names: List[str]
    indexs: Tuple[int, int]
    sequence: str

# class Vertice:
#     name:
#     chain:
#


def extract_sequence_with_seqio(mmcif_path,af_version: int):
    """
    Extracts the sequence from an mmCIF file using Bio.SeqIO.

    Args:
        mmcif_path (str): Path to the mmCIF file.
        af_version if 2 then PDB and if 3 cif

    Returns:
        str: The amino acid sequence as a single-letter code string.
    """
    format = {'2':"pdb-atom", '3':"cif-atom"}
    sequences = []
    for record in SeqIO.parse(mmcif_path, format[af_version]):
        sequences.append(str(record.seq))
        print(record.id)
    return ''.join(sequences)


def extract_sequence_from_mmcif(mmcif_path):
    """
    Extracts the amino acid sequence from the ATOM records of an mmCIF file.

    Args:
        mmcif_path (str): Path to the mmCIF file.

    Returns:
        str: The amino acid sequence as a single-letter code string.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("Model", mmcif_path)

    # Use the PPBuilder to build polypeptides and extract the sequence
    ppb = PPBuilder()
    sequences = []

    for pp in ppb.build_peptides(structure):
        sequences.append(pp.get_sequence())  # Seq object
        pp.get_id()

    # Combine sequences if there are multiple chains
    return ''.join(str(seq) for seq in sequences)


def find_high_confidence_regions(plddt_array, confidence_threshold=40, gap_threshold=1):
    """
    Finds ranges of high-confidence regions in a plDDT array while preserving the order.

    Args:
        plddt_array (np.ndarray): Array containing plDDT values for each residue.
        confidence_threshold (float): Minimum confidence value to include a residue.
        gap_threshold (int): Maximum allowed gap between indices in the same region.

    Returns:
        list[tuple[int, int]]: A list of tuples, each representing the start and end indices of a high-confidence region.
    """
    # Find indices where plDDT is above the confidence threshold
    indices = np.where(plddt_array > confidence_threshold)[0]

    if len(indices) == 0:
        return []  # No high-confidence regions

    regions = []
    start_index = indices[0]

    for i in range(1, len(indices)):
        if indices[i] - indices[i - 1] > gap_threshold:
            # Add the current region when a gap is found
            regions.append((int(start_index), int(indices[i - 1])))
            start_index = indices[i]

    # Append the last region
    regions.append((int(start_index), int(indices[-1])))
    return regions


def extract_subunit_info(indexs: List[Tuple[int, int]], token_chain_ids: List[str], full_seq: str) -> List[SubunitInfo]:
    subunit_infos = []
    # todo: change subunit name to be from file
    chain_occ_counter = Counter()  # Counter to track occurrences of each chain ID

    for start, end in indexs:
        # Extract chain IDs in the current node and ensure uniqueness
        chains_ids_in_node = list(dict.fromkeys(token_chain_ids[start:end + 1]))  # keep order
        subunit_name = "".join(f"{chain_id}{chain_occ_counter[chain_id] + 1}" for chain_id in chains_ids_in_node)
        for chain_id in chains_ids_in_node:
            chain_occ_counter[chain_id] += 1
        subunit_infos.append(SubunitInfo(
            name=subunit_name,
            chain_names=chains_ids_in_node,
            indexs=(start, end),
            sequence=full_seq[start:end + 1]
        ))

    return subunit_infos


def atom_plddt_to_res_plddt(structure, atom_plddts):
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

# for each multimer in Ben's data print how many nodes, len before and after and


def plot_pae_plddt(pae_as_arr: np.array, plddt_array, nodes, edges, plot_name: str): #todo: add edge width due weight
    # Define AlphaFold plDDT color scheme
    plddt_cmap = ListedColormap(['#FF7D45', '#FFDB13', '#65CBF3', '#0053D6'])  # Orange, Yellow, Cyan, Blue
    bounds = [0, 50, 70, 90, 100]  # Define confidence levels
    norm = BoundaryNorm(bounds, plddt_cmap.N)  # Use BoundaryNorm for discrete segments

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(15, 15))
    matrix_ax = ax.matshow(pae_as_arr, vmin=0., vmax=np.max(pae_as_arr), cmap='Greens_r')
    fig.colorbar(matrix_ax, label='PAE', ax=ax)

    # Add plDDT bars to X and Y axes
    divider_width = 0.02
    divider_offset = 0.02

    # plDDT bar for X-axis
    x_cb_ax = fig.add_axes([
        ax.get_position().x0,
        ax.get_position().y1 + divider_offset,  # Position above the matrix
        ax.get_position().width,
        divider_width
    ])
    x_cb_ax.imshow(
        plddt_array.reshape(1, -1),
        aspect='auto',
        cmap=plddt_cmap,
        norm=norm
    )
    x_cb_ax.set_xticks([])
    x_cb_ax.set_yticks([])

    # plDDT bar for Y-axis
    y_cb_ax = fig.add_axes([
        ax.get_position().x1 + divider_offset,  # Position to the right of the matrix
        ax.get_position().y0,
        divider_width,
        ax.get_position().height
    ])
    y_cb_ax.imshow(
        plddt_array.reshape(-1, 1),
        aspect='auto',
        cmap=plddt_cmap,
        norm=norm
    )
    y_cb_ax.set_xticks([])
    y_cb_ax.set_yticks([])

    # Draw square borders and annotate with names
    for name, start, end in nodes:
        rect = Rectangle(
            (start - 0.5, start - 0.5),  # Bottom-left corner
            end - start,  # Width
            end - start,  # Height
            edgecolor='blue',
            facecolor='none',
            linewidth=2
        )
        ax.add_patch(rect)
        # Annotate the square
        ax.text(
            start + (end - start) / 2,  # X position (center of the square)
            start + (end - start) / 2,  # Y position (center of the square)
            name,  # Annotation text
            color='black',
            fontsize=10,
            ha='center',
            va='center',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.6)  # Optional: Background for text
        )

    # from here code for edges
    # Mapping of square names to their centers for drawing edges
    square_centers = {
        name: (start + (end - start) / 2, start + (end - start) / 2)  # Center of the square
        for name, start, end in nodes
    }

    # Draw edges as L-shaped lines
    for edge in edges:
        if edge[0] in square_centers and edge[1] in square_centers:
            x1, y1 = square_centers[edge[0]]
            x2, y2 = square_centers[edge[1]]
            # Break the L-shape into two parts
            # Horizontal line
            ax.plot([x1, x2], [y1, y1], color='#800080', linestyle='-', linewidth=2.5)
            # Vertical line
            ax.plot([x2, x2], [y1, y2], color='#800080', linestyle='-', linewidth=2.5)

    # Add a plDDT legend
    legend_elements = [
        Patch(facecolor='#FF7D45', edgecolor='black', label='Very low (0-50)'),
        Patch(facecolor='#FFDB13', edgecolor='black', label='Low (50-70)'),
        Patch(facecolor='#65CBF3', edgecolor='black', label='Confident (70-90)'),
        Patch(facecolor='#0053D6', edgecolor='black', label='Very high (90-100)')
    ]
    ax.legend(
        handles=legend_elements,
        loc='best',  # Position the legend above the plot
        bbox_to_anchor=(1, 1.2),  # Center above the plot with some padding
        title="plDDT Confidence",
        frameon=True
    )

    # Add title above the entire plot
    fig.suptitle("Predicted Aligned Error (PAE) with plDDT Bars \nand Confidence>40 Borders", fontsize=16, x=0.4, y=0.9)

    # Show the plot
    # plt.savefig('all_plot.png', format='png', dpi=300, bbox_inches='tight')
    plt.savefig(plot_name + 'all_plot.png', format='png', dpi=300, bbox_inches='tight')
    plt.show()


def plot_pae_plddt2(pae_as_arr: np.array, plddt_array, nodes, edges, plot_name: str, use_edge_weights: bool = True, weight_threshold: float = 15):
    """
    Plot PAE matrix with plDDT confidence bars and optional edge thickness based on weights.

    Args:
        pae_as_arr (np.array): Predicted Aligned Error (PAE) matrix.
        plddt_array (np.array): Array of plDDT confidence scores.
        nodes (list[tuple]): List of tuples representing nodes (name, start, end).
        edges (list[tuple]): List of tuples representing edges (node1, node2, weight).
        plot_name (str): Name of the file to save the plot.
        use_edge_weights (bool): If True, edge thickness reflects weights.
        weight_threshold (float): Maximum weight value for normalization.
    """
    # Define AlphaFold plDDT color scheme
    plddt_cmap = ListedColormap(['#FF7D45', '#FFDB13', '#65CBF3', '#0053D6'])
    bounds = [0, 50, 70, 90, 100]
    norm = BoundaryNorm(bounds, plddt_cmap.N)

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(15, 15))
    matrix_ax = ax.matshow(pae_as_arr, vmin=0., vmax=np.max(pae_as_arr), cmap='Greens_r')
    fig.colorbar(matrix_ax, label='PAE', ax=ax)

    # Add plDDT bars to X and Y axes
    divider_width = 0.02
    divider_offset = 0.02

    # X-axis plDDT bar
    x_cb_ax = fig.add_axes([
        ax.get_position().x0,
        ax.get_position().y1 + divider_offset,
        ax.get_position().width,
        divider_width
    ])
    x_cb_ax.imshow(plddt_array.reshape(1, -1), aspect='auto', cmap=plddt_cmap, norm=norm)
    x_cb_ax.set_xticks([])
    x_cb_ax.set_yticks([])

    # Y-axis plDDT bar
    y_cb_ax = fig.add_axes([
        ax.get_position().x1 + divider_offset,
        ax.get_position().y0,
        divider_width,
        ax.get_position().height
    ])
    y_cb_ax.imshow(plddt_array.reshape(-1, 1), aspect='auto', cmap=plddt_cmap, norm=norm)
    y_cb_ax.set_xticks([])
    y_cb_ax.set_yticks([])

    # Draw square borders and annotate nodes
    square_centers = {}
    for name, start, end in nodes:
        rect = Rectangle((start - 0.5, start - 0.5), end - start, end - start, edgecolor='blue', facecolor='none', linewidth=2)
        ax.add_patch(rect)
        ax.text(start + (end - start) / 2, start + (end - start) / 2, name, color='black', fontsize=10, ha='center', va='center',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.6))
        square_centers[name] = (start + (end - start) / 2, start + (end - start) / 2)

    # Normalize edge weights
    if use_edge_weights:
        max_thickness = 5.0
        min_thickness = 0.5
        normalized_weights = [
            (edge[0], edge[1], max(min((edge[2] / weight_threshold) * (max_thickness - min_thickness) + min_thickness, max_thickness), min_thickness))
            for edge in edges
        ]
        print(normalized_weights)
    else:
        normalized_weights = [(edge[0], edge[1], 2.0) for edge in edges]  # Default thickness

    # Draw edges
    for edge in normalized_weights:
        if edge[0] in square_centers and edge[1] in square_centers:
            x1, y1 = square_centers[edge[0]]
            x2, y2 = square_centers[edge[1]]
            ax.plot([x1, x2], [y1, y1], color='#800080', linestyle='-', linewidth=edge[2])  # Horizontal line
            ax.plot([x2, x2], [y1, y2], color='#800080', linestyle='-', linewidth=edge[2])  # Vertical line

    # Add plDDT legend
    legend_elements = [
        Patch(facecolor='#FF7D45', edgecolor='black', label='Very low (0-50)'),
        Patch(facecolor='#FFDB13', edgecolor='black', label='Low (50-70)'),
        Patch(facecolor='#65CBF3', edgecolor='black', label='Confident (70-90)'),
        Patch(facecolor='#0053D6', edgecolor='black', label='Very high (90-100)')
    ]
    ax.legend(handles=legend_elements, loc='best', bbox_to_anchor=(1, 1.2), title="plDDT Confidence", frameon=True)

    # Add title and save
    fig.suptitle("Predicted Aligned Error (PAE) with plDDT Bars and Edge Weights", fontsize=16, x=0.4, y=0.9)
    plt.savefig(f'{plot_name}_plot.png', format='png', dpi=300, bbox_inches='tight')
    plt.show()


def find_edges(subunits_info: SubunitInfo, pae_matrix: np.array, threshold: int = 15) -> list[tuple[str, str, float]]:
    edges = []
    for subunit1, subunit2 in itertools.combinations(subunits_info, 2):
        pae_rect = pae_matrix[subunit1.indexs[0]:subunit1.indexs[1], subunit2.indexs[0]:subunit2.indexs[1]]
        # pae_rect.
        # print(f"[{node1[0]}:{node1[1]][node2[0]:node2[1]]")
        pae_score = np.mean(pae_rect)
        # print(pae_score)
        if pae_score < threshold:
            edges.append((subunit1.name, subunit2.name, float(pae_score)))
    return edges


def get_chain_ids_per_residue(structure):
    """
    Made for getting the token_chain_ids which require for extract_subunit_info() in AF2 case. #todo: not a good practice
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

def graph(structure_path: str, data_path:str, af_version: str)->tuple[list,list]:
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
        plddt_array = np.array(json_full_data['pae'])  # todo: not sure if it keeps order this way
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("original_pdb", structure_path)
        token_chain_ids = get_chain_ids_per_residue(structure)
    full_seq = extract_sequence_with_seqio(structure_path,
                                           af_version)  # todo: take the full seq from complex file instead!!
    groups_indexs = find_high_confidence_regions(plddt_array, confidence_threshold=40)

    filename = os.path.basename(data_path)  # Extract 'cdf_ddf' #todo: huraney
    parts = filename.split('_')  # Now split by '_'
    replacement_dict = {'A': parts[0], 'B': parts[1]}  # Use only the last part
    token_chain_ids_updated = [replacement_dict.get(item, item) for item in token_chain_ids]

    subunits_info = extract_subunit_info(groups_indexs, token_chain_ids_updated, full_seq)
    vertices = []
    for subunit in subunits_info:
        vertices.append({'name': subunit.name, 'chain': subunit.chain_names[0],
                         'start': subunit.indexs[0], 'end': subunit.indexs[1]})
    edges = find_edges(subunits_info, pae_as_arr, threshold=15)
    return (vertices, edges)

if __name__ == '__main__':
    if len(sys.argv) == 4:
        structure_path, data_path, af_version = os.path.abspath(sys.argv[1]),os.path.abspath(sys.argv[2]),sys.argv[3]
    else:
        print("usage: <script> structure_path data_path af_version")
    g = graph(structure_path, data_path, af_version)
    print (g)
    #plot_pae_plddt(pae_as_arr, plddt_array, nodes_as_req, edges, 'skip4_pae15_')
    #plot_pae_plddt2(pae_as_arr, plddt_array, nodes_as_req, edges, 'with_weights')
