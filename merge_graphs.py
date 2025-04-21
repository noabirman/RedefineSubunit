import copy
import json
import os
import sys
from dataclasses import asdict

import networkx as nx
from complex_graph import graph, SubunitInfo, extract_sequence_with_seqio
from typing import List


def overlap(v1: SubunitInfo, v2: SubunitInfo) -> bool:
    """Check if two SubunitInfo nodes overlap in at least one chain."""
    return any(chain in v2.chain_names for chain in v1.chain_names) and not (v1.end < v2.start or v2.end < v1.start)


def merge_graphs(graphs: List[nx.Graph], folder: str) -> nx.Graph:
    """Merge nodes with overlapping indices in the same chain across multiple graphs."""
    # Build an overlap graph
    overlap_graph = create_overlap_graph(graphs)
    # Merge nodes in each component
    merged_graph  = merge_connected_components(overlap_graph)


    # save subunits to json
    with open(os.path.join(os.path.dirname(folder), "subunits.json"), 'w') as f:
        json.dump(output_subunits, f, indent=2)
    # Create a deep copy of the merged graph
    merged_graph_copy = copy.deepcopy(merged_graph)

    # Convert SubunitInfo objects to dictionaries in the copied graph
    for node in merged_graph_copy.nodes(data=True):
        if 'data' in node[1] and isinstance(node[1]['data'], SubunitInfo):
            node[1]['data'] = asdict(node[1]['data'])
    # create json file
    with open(os.path.join(os.path.dirname(folder), "merged_graph.json"), 'w') as f:
        json.dump(nx.readwrite.json_graph.node_link_data(merged_graph_copy), f, indent=2)

    return merged_graph


def create_graphs_from_folder(folder_path):
    print(f"Creating graphs from folder {folder_path}")
    graphs = []
    # Iterate over all folders in the folder
    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)
        if os.path.isdir(item_path):
            print(f"Processing {item}")
            data_path = os.path.join(item_path, f"{item}_confidences.json")
            data_path = os.path.abspath(data_path)
            structure_path = os.path.join(item_path, f"{item}_model.cif")
            structure_path = os.path.abspath(structure_path)
            if not os.path.exists(data_path) or not os.path.exists(structure_path):
                print(f"Skipping {item}: Required files not found.")
                continue
            graph1 = graph(structure_path, data_path, af_version='3')
            graphs.append(graph1)
    return graphs


def create_overlap_graph(graphs):
    overlap_graph = nx.Graph()  # Temporary graph for finding connected components

    # Add all nodes to the overlap graph
    for idx, G in enumerate(graphs):
        for node, data in G.nodes(data=True):
            unique_node = f"{node}_{idx}"
            overlap_graph.add_node(unique_node, data=data['data'])

    # Add edges between overlapping nodes
    nodes_list = list(overlap_graph.nodes.items())  # Convert to list for pairwise comparison
    for i in range(len(nodes_list)):
        name1, data1 = nodes_list[i]
        subunit1 = data1['data']
        for j in range(i + 1, len(nodes_list)):
            name2, data2 = nodes_list[j]
            subunit2 = data2['data']
            if overlap(subunit1, subunit2):
                overlap_graph.add_edge(name1, name2)

    return overlap_graph


def merge_connected_components(overlap_graph):
    connected_components = list(nx.connected_components(overlap_graph))
    merged_graph = nx.Graph()
    node_mapping = {}
    node_dict = {name: data['data'] for name, data in overlap_graph.nodes(data=True)}
    for component in connected_components:
        subunits = [node_dict[member] for member in component]
        # Merge properties
        merged_name = f"{list(component)[0]}_high"  # Use the first node's name as the merged name
        merged_chains = sorted(set(chain for subunit in subunits for chain in subunit.chain_names))  # change
        merged_start = min(subunit.start for subunit in subunits)
        merged_end = max(subunit.end for subunit in subunits)
        merged_sequence = merge_names(subunits, merged_start, merged_end)
        merged_subunit = SubunitInfo(name=merged_name, chain_names=merged_chains, start=merged_start, end=merged_end,
                                     sequence=merged_sequence)
        # Add merged node to the merged graph
        merged_graph.add_node(merged_name, data=merged_subunit)
        # output_subunits[merged_name] = {
        #     "name": merged_name,
        #     "chain_names": merged_chains,
        #     "start_res": merged_start,
        #     "sequence": merged_sequence,
        # }
        # Map original nodes to merged node name
        for name in component:
            node_mapping[name] = merged_name
        for idx, G in enumerate(graphs):
            for u, v in G.edges():
                merged_u = node_mapping[f"{u}_{idx}"]
                merged_v = node_mapping[f"{v}_{idx}"]
                merged_graph.add_edge(merged_u, merged_v)
    return merged_graph


def merge_names(subunits, start, end):
    merged_sequence = [''] * (end + 1 - start)  # inclusion of end position
    for subunit in subunits:
        for i, char in enumerate(subunit.sequence):
            pos = subunit.start - start + i
            if merged_sequence[pos] == '':
                merged_sequence[pos] = char
            elif merged_sequence[pos] != char:
                raise ValueError(f"Conflict detected at position {pos}: {merged_sequence[pos]} vs {char}")
    merged_sequence = "".join(merged_sequence)
    return merged_sequence

def save_subunits_info(graph, original_subunits_path):

# main
if __name__ == "__main__":
    if len(sys.argv) == 3:
        folder_path = os.path.abspath(sys.argv[1])
        original_subunits_path = os.path.abspath(sys.argv[2])
        graphs = create_graphs_from_folder(folder_path)
        merged_graph = merge_graphs(graphs, folder_path)
        save_subunits_info(merged_graph, original_subunits_path)
        # nodes_with_edges = [node for node, degree in merged_graph.degree() if degree > 0]
        # sub_graph = merged_graph.subgraph(nodes_with_edges).copy()

    else:
        print("usage: <script> enter folder_name, subunits_json_path")
