import json
import os
import sys
from collections import defaultdict

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
    merged_graph  = merge_connected_components(overlap_graph, graphs)

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

def merge_connected_components(overlap_graph, graphs: List[nx.Graph]):
    connected_components = list(nx.connected_components(overlap_graph))
    merged_graph = nx.Graph()
    node_mapping = {}
    node_dict = {name: data['data'] for name, data in overlap_graph.nodes(data=True)}
    for index, component in enumerate( connected_components):
        subunits = [node_dict[member] for member in component]
        # Merge properties
        merged_name = f"{list(component)[0]}_high"  # Use the first node's name as the merged name
        merged_chains = sorted(set(chain for subunit in subunits for chain in subunit.chain_names))  # change
        merged_start = min(subunit.start for subunit in subunits)
        merged_end = max(subunit.end for subunit in subunits)
        merged_sequence = merge_sequence(subunits, merged_start, merged_end)
        merged_subunit = SubunitInfo(name=merged_name, chain_names=merged_chains, start=merged_start, end=merged_end,
                                     sequence=merged_sequence)
        # Add merged node to the merged graph
        merged_graph.add_node(merged_name, data=merged_subunit)
        # Map original nodes to merged node name
        for name in component:
            node_mapping[name] = merged_name
    for idx, G in enumerate(graphs):
        for u, v in G.edges():
            merged_u = node_mapping[f"{u}_{idx}"]
            merged_v = node_mapping[f"{v}_{idx}"]
            merged_graph.add_edge(merged_u, merged_v)
    return merged_graph

def merge_sequence(subunits, start, end):
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

def save_subunits_info(graph: nx.Graph, subunit_name_mapping_path: str, subunits_info_path: str, folder) -> None:
    """
    Process graph nodes (high segments) and create a unified JSON with both high and low segments.
    Node data in graph should be SubunitInfo objects.
    """
    # Load subunit name mapping
    with open(subunit_name_mapping_path, 'r') as f:
        name_mapping = json.load(f)

    # Load subunits info
    with open(subunits_info_path, 'r') as f:
        subunits_info = json.load(f)

    # Dictionary to store all segments (high and low)
    unified_subunits = {}

    # Group high segments by original subunit
    high_segments_by_subunit = defaultdict(list)

    # Process high segments from graph
    for node in graph.nodes:
        # Extract base name from node name (e.g., "A" from "A1_high")
        base_name = node.split('_')[0]
        # Get original subunit name from mapping
        original_chain_name  = name_mapping[base_name]
        # Find the original_name where subunits_info[key]['chain_names'] contains original_chain_name
        original_name = next(
            (key for key, info in subunits_info.items() if original_chain_name in info['chain_names']),
            None
        )

        if original_name is None:
            raise ValueError(f"No subunit found with chain name: {original_chain_name}")
        # Find the subunit_info where one of the chain_names matches the original_name
        subunit_info = subunits_info.get(original_name)

        if subunit_info is None:
            raise ValueError(f"No subunit found with chain name: {original_name}")
        # added for debug
        print("→ original_name chosen:", original_name)
        print("→ JSON chain_names:", subunit_info['chain_names'])
        print("→ JSON sequence length:", len(subunit_info['sequence']))
        print("→ JSON sequence start:", subunit_info['sequence'][:20], "…")
        # Get node data (SubunitInfo object)
        node_data = graph.nodes[node]['data']
        start = node_data.start
        end = node_data.end  # Remember: this is inclusive
        sequence = node_data.sequence

        # Verify sequence matches the original
        #expected_sequence = subunit_info['sequence'][start-1:end]
        # Extract what we think is the “correct” slice
        full_seq = subunit_info['sequence']
        # If your SubunitInfo.start/end are *inclusive* 1-based, this is right:
        expected_sequence = full_seq[start - 1:end]
        # But if end were already exclusive, you’d want full_seq[start-1:end-1]…

        if sequence != expected_sequence:
            print("==== DEBUG SEQUENCE MISMATCH ====")
            print(f"Node:           {node}")
            print(f"Start, end:     {start}, {end}")
            print(f"SubunitInfo.seq:   {sequence!r}  (len={len(sequence)})")
            print(f"Expected slice:    {expected_sequence!r}  (len={len(expected_sequence)})")
            print(f"Full sequence len: {len(full_seq)}")
            print(f"First 10 res:      {full_seq[:10]!r}")
            print(f"Chain names:       {subunit_info['chain_names']!r}")
            print(f"Name mapping JSON: {name_mapping}")
            print(f"Subunits_info keys:{list(subunits_info.keys())[:10]} …")
            raise ValueError("Stopping here for debug")
        # if sequence != expected_sequence:
        #     raise ValueError(f"Sequence mismatch in node '{node}'. Expected: {expected_sequence}, Found: {sequence}")

        # Store high segment information
        high_segments_by_subunit[original_name].append({
            'node_name': node,
            'start': start,
            'end': end
        })

        # Create high segment entry
        unified_subunits[node] = {
            'name': node,
            'sequence': sequence,
            'chain_names': subunit_info['chain_names'],
            'start_res': start + 1  # Convert to 1-based indexing
        }

    # Process each original subunit to create low segments
    for original_name, segments in high_segments_by_subunit.items():
        subunit_info = subunits_info[original_name]
        full_sequence = subunit_info['sequence']
        total_length = len(full_sequence)

        # Sort segments by start position
        segments.sort(key=lambda x: x['start'])

        # Track the end of the last processed segment
        last_end = -1
        low_segment_index = 1

        # Get base name from any segment's node name
        base_name = segments[0]['node_name'].split('_')[0]
        # Process gaps before first high segment and between high segments
        for segment in segments:
            current_start = segment['start']

            # If there's a gap before this segment
            if current_start > last_end + 1:
                low_start = last_end + 1
                low_end = current_start - 1

                # Create low segment name
                low_name = f"{base_name}_low_{low_segment_index}"

                # Extract sequence for low segment
                low_sequence = full_sequence[low_start:low_end + 1]  # +1 because end is inclusive

                # Create low segment entry
                unified_subunits[low_name] = {
                    'name': low_name,
                    'sequence': low_sequence,
                    'chain_names': subunit_info['chain_names'],
                    'start_res': low_start + 1  # Convert to 1-based indexing
                }

                low_segment_index += 1

            last_end = segment['end']

        # Check for gap after last high segment
        if last_end < total_length - 1:
            low_start = last_end + 1
            low_end = total_length - 1

            # Create low segment name
            low_name = f"{base_name}_low_{low_segment_index}"

            # Extract sequence for low segment
            low_sequence = full_sequence[low_start:low_end + 1]  # +1 because end is inclusive

            # Create low segment entry
            unified_subunits[low_name] = {
                'name': low_name,
                'sequence': low_sequence,
                'chain_names': subunit_info['chain_names'],
                'start_res': low_start + 1  # Convert to 1-based indexing
            }

    sorted_unified_subunits = dict(
        sorted(
            unified_subunits.items(),
            key=lambda item: (item[1]['name'].split('_')[0], item[1]['start_res'])
        )
    )

    output_folder = os.path.join(os.path.dirname(folder), 'combfold')
    os.makedirs(output_folder, exist_ok=True)
    output_json_path = os.path.join(output_folder, 'subunits_info.json')

    with open(output_json_path, 'w') as f:
        json.dump(sorted_unified_subunits, f, indent=4)
# main
if __name__ == "__main__":
    if len(sys.argv) >= 2:
        folder_path = os.path.abspath(sys.argv[1])
        original_subunits_path = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.join(os.path.dirname(folder_path), 'subunits_info.json')
        mapping_path = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.path.join(os.path.dirname(folder_path), 'chain_id_mapping.json')

        graphs = create_graphs_from_folder(folder_path)
        merged_graph = merge_graphs(graphs, folder_path)
        save_subunits_info(merged_graph, mapping_path, original_subunits_path, folder_path)
    else:
        print("usage: <script> enter folder_name, [subunits_json_path], [mapping_json_path], [output_json_path]")
        #