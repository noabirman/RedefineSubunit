import json
import os
import sys
from collections import defaultdict
import difflib
import networkx as nx
from create_graph_from_af_model import graph, SubunitInfo
from typing import List

from results_analysis.vizualization_plots import show_circle, show_subunit_circle_graph
import pickle


def check_subunit_sequence_reconstruction(original_path, new_path):
    import json

    def load(path):
        with open(path) as f:
            return json.load(f)

    original_data = load(original_path)
    new_data = load(new_path)

    for orig_id, orig_info in original_data.items():
        matching_subs = [
            (k, v) for k, v in new_data.items()
            if k.startswith(orig_id + "_")
        ]

        if not matching_subs:
            raise ValueError(f"[Missing] No matching sub-subunits for original subunit '{orig_id}'.")

        # Sort sub-subunits by start_res
        matching_subs.sort(key=lambda kv: kv[1]["start_res"])

        expected_start = orig_info["start_res"]
        expected_chains = set(orig_info["chain_names"])

        reconstructed_seq = ""

        for k, v in matching_subs:
            start = v["start_res"]
            seq = v["sequence"]
            end = start + len(seq) - 1

            if start != expected_start:
                raise ValueError(
                    f"[Gap/Overlap] Subunit '{k}' expected to start at residue {expected_start}, but starts at {start}."
                )

            if start == 1 and (not seq or seq[0] != "M"):
                raise ValueError(
                    f"[BadStart] Subunit '{k}' starts at residue 1 but first amino acid is not 'M': '{seq[:10]}'"
                )

            if set(v["chain_names"]) != expected_chains:
                raise ValueError(
                    f"[Chain mismatch] Subunit '{k}' has chains {v['chain_names']}, expected {orig_info['chain_names']}."
                )

            reconstructed_seq += seq
            expected_start = end + 1

        if reconstructed_seq != orig_info["sequence"]:
            diff = '\n'.join(difflib.ndiff(orig_info["sequence"], reconstructed_seq))
            raise ValueError(
                f"[Mismatch] Reconstructed sequence for '{orig_id}' does not match original.\n"
                f"  Original     : {orig_info['sequence']}\n"
                f"  Reconstructed: {reconstructed_seq}"
                f"  Diff:\n{diff}"
            )

def overlap(v1: SubunitInfo, v2: SubunitInfo, threshold = 5) -> bool:
    """Check if two SubunitInfo nodes overlap in at least one chain."""
    return any(chain in v2.chain_names for chain in v1.chain_names) and not (v1.end + threshold < v2.start or v2.end + threshold < v1.start)

def merge_graphs(graphs: List[nx.Graph], name_mapping: object, subunits_info: object) -> nx.Graph:
    """Merge nodes with overlapping indices in the same chain across multiple graphs."""
    # Build an overlap graph
    overlap_graph = create_overlap_graph(graphs)
    # Merge nodes in each component
    merged_graph  = merge_connected_components(overlap_graph, graphs,name_mapping,subunits_info)

    return merged_graph

# def create_graphs_from_folder(folder_path):
#     print(f"Creating graphs from folder {folder_path}")
#     graphs = []
#     # Iterate over all folders in the folder
#     for item in os.listdir(folder_path):
#         item_path = os.path.join(folder_path, item)
#         if os.path.isdir(item_path):
#             print(f"Processing {item}") # todo: print this I just removed it for debug
#             data_path = os.path.join(item_path, f"{item}_confidences.json")
#             data_path = os.path.abspath(data_path)
#             structure_path = os.path.join(item_path, f"{item}_model.cif")
#             structure_path = os.path.abspath(structure_path)
#             if not os.path.exists(data_path) or not os.path.exists(structure_path):
#                 print(f"Skipping {item}: Required files not found.")
#                 continue
#             graph1 = graph(structure_path, data_path, af_version='3')
#             graphs.append(graph1)
#     return graphs
def create_graphs_from_folder(folder_path, use_ipsae=False):
    print(f"Creating graphs from folder {folder_path} (Mode: {'ipSAE' if use_ipsae else 'pLDDT'})")
    graphs = []

    # Define ipsae folder path (sibling directory)
    ipsae_results_dir = os.path.join(os.path.dirname(folder_path), "ipsae_results")

    # Cutoffs must match your ipsae run
    pae_cut = "15"
    dist_cut = "10"

    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)
        if os.path.isdir(item_path):
            print(f"Processing {item}")
            data_path = os.path.join(item_path, f"{item}_confidences.json")
            structure_path = os.path.join(item_path, f"{item}_model.cif")

            # Absolute paths
            data_path = os.path.abspath(data_path)
            structure_path = os.path.abspath(structure_path)

            if not os.path.exists(data_path) or not os.path.exists(structure_path):
                print(f"Skipping {item}: Required files not found.")
                continue

            ipsae_path = None
            if use_ipsae:
                # 1. Check if the main results folder exists
                if not os.path.exists(ipsae_results_dir):
                    raise FileNotFoundError(
                        f"Error: You requested ipSAE mode, but the folder '{ipsae_results_dir}' does not exist.")

                # 2. Find the specific file
                expected_name = f"{item}_{item}_model_{pae_cut}_{dist_cut}_byres.txt"
                candidate_path = os.path.join(ipsae_results_dir, expected_name)

                # Fallback search
                if not os.path.exists(candidate_path):
                    candidates = [f for f in os.listdir(ipsae_results_dir) if
                                  f.startswith(item) and f.endswith("_byres.txt")]
                    if candidates:
                        candidate_path = os.path.join(ipsae_results_dir, candidates[0])

                if os.path.exists(candidate_path):
                    ipsae_path = candidate_path
                else:
                    # STRICT BEHAVIOR: Raise error if missing, so you know the comparison failed
                    print(f"Warning: No ipSAE file found for {item}. Skipping this graph.")
                    continue

                    # Pass the path (or None) to the graph function
            graph1 = graph(structure_path, data_path, af_version='3', ipsae_path=ipsae_path)
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

def merge_connected_components(overlap_graph, graphs: List[nx.Graph],subunit_name_mapping_path: str, subunits_info_path: str):
    connected_components = list(nx.connected_components(overlap_graph))
    merged_graph = nx.Graph()
    node_mapping = {}
    node_dict = {name: data['data'] for name, data in overlap_graph.nodes(data=True)}
    # for name, sub in node_dict.items():
    #     print("DEBUG node:", name, sub.start, sub.end, sub.chain_names)

    # Group components by chain prefix
    chain_components = defaultdict(list)
    for component in connected_components:
        first_node = next(iter(component))
        chain_prefix = first_node.split('_')[0]

        # Verify all nodes have same prefix
        if not all(node.split('_')[0] == chain_prefix for node in component):
            raise ValueError(f"Component contains nodes with different chain prefixes: {component}")

        chain_components[chain_prefix].append(component)

    # Sort components within each chain by start position
    for chain_prefix in chain_components:
        chain_components[chain_prefix].sort(
            key=lambda comp: min(node_dict[member].start for member in comp)
        )

    # Process each chain's sorted components
    for chain_prefix, components in chain_components.items():
        for index, component in enumerate(components):
            is_last = index == len(components) - 1  # True if this is the last component
            subunits = [node_dict[member] for member in component]
            # Merge properties
            merged_name = f"{chain_prefix}_high_{index+1}"  # Use the first node's name as the merged name
            merged_chains = sorted(set(chain for subunit in subunits for chain in subunit.chain_names))  # change
            merged_start = min(subunit.start for subunit in subunits)
            if merged_start <= 5:
                merged_start = 1
            merged_end = max(subunit.end for subunit in subunits)
            if is_last:
                _, subunit = find_original_subunit_info(chain_prefix, name_mapping, subunits_info, merged_start, merged_end)
                if merged_end > (len(subunit["sequence"]) - 5):
                    merged_end = len(subunit["sequence"])
            merged_sequence = merge_sequence(subunits, merged_start, merged_end)
            merged_subunit = SubunitInfo(name=merged_name, chain_names=merged_chains, start=merged_start, end=merged_end,
                                         sequence=merged_sequence)
            print(f"→ Merged node: {merged_name}, chains={merged_chains}, start={merged_start}, end={merged_end}, seq_len={len(merged_sequence)}")

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
    merged_sequence = ['-'] * (end + 1 - start)  # inclusion of end position
    print(f"→ Merging sequence from {start} to {end}")
    for subunit in subunits:
        print(f"  ↳ Subunit: {subunit.name}, start={subunit.start}, end={subunit.end}, seq_len={len(subunit.sequence)} seq={subunit.sequence}")
        for i, char in enumerate(subunit.sequence):
            pos = subunit.start - start + i
            if pos < 0 or pos >= len(merged_sequence):
                raise IndexError(
                    f"Index out of bounds while merging: pos={pos}, start={start}, subunit.start={subunit.start}, i={i}, merged_len={len(merged_sequence)}")
            if merged_sequence[pos] == '-':
                merged_sequence[pos] = char
            elif merged_sequence[pos] != char:
                print(f"Conflict at pos={pos + start}: existing={merged_sequence[pos]} vs new={char}")
                raise ValueError(f"Conflict detected at position {pos}: {merged_sequence[pos]} vs {char}")
    merged_sequence = "".join(merged_sequence)
    return merged_sequence

def sequences_match(seq1: str, seq2: str) -> bool:
    return all(a == b or a == '-' or b == '-' for a, b in zip(seq1, seq2)) and len(seq1) == len(seq2)

def find_original_subunit_info(base_name: str, name_mapping: dict, subunits_info: dict, start, end) -> tuple[str, dict]:
    """
    Find the original subunit name and info given a base name using chain name mapping.

    Parameters
    ----------
    base_name : str
        The base name to look up (e.g., "A" from "A_high_1")
    name_mapping : dict
        Mapping of base names to chain names
    subunits_info : dict
        Dictionary containing subunit information with chain names

    Returns
    -------
    tuple[str, dict]
        Original subunit name and its info dictionary

    Raises
    ------
    ValueError
        If no subunit is found with the mapped chain name
    """
    original_chain_name = name_mapping[base_name]['chain_id']
    matching_subunits = {
        key: info for key, info in subunits_info.items()
        if original_chain_name in info['chain_names']
    }

    if not matching_subunits:
        raise ValueError(f"No subunit found with chain name: {original_chain_name}")

    if len(matching_subunits) == 1:
        # If only one match, return it directly
        original_subunit_name = next(iter(matching_subunits))
    else:
        # If multiple matches, filter by start/end position
        raise ValueError(f"Multiple subunits found with chain name: {original_chain_name}")

    if original_subunit_name is None:
        raise ValueError(f"No subunit found matching position criteria for chain: {original_chain_name}")

    subunit_info = subunits_info[original_subunit_name]
    return original_subunit_name, subunit_info

def save_subunits_info(graph: nx.Graph, name_mapping: dict, subunits_info: dict, folder: str,
                           output_subdir: str = 'combfold') -> None:
    """
    Process graph nodes (high segments) and create a unified JSON with both high and low segments.

    Parameters
    ----------
    graph : nx.Graph
        Graph with SubunitInfo objects as node data
    name_mapping : dict
        Mapping of base names to chain names
    subunits_info : dict
        Dictionary containing subunit information
    folder : str
        Path to output folder
    output_subdir: str
    """
    output_folder = os.path.join(os.path.dirname(folder), output_subdir)

    # Dictionary to store all segments (high and low)
    unified_subunits = {}

    # Group high segments by original subunit
    high_segments_by_subunit = defaultdict(list)

    # Process high segments from graph
    for node in graph.nodes:
        # Extract base name from node name (e.g., "A" from "A1_high")
        base_name = node.split('_')[0]
        # Get original subunit name from mapping
        node_data = graph.nodes[node]['data']
        original_name, subunit_info = find_original_subunit_info(base_name, name_mapping, subunits_info, node_data.start, node_data.end)
        # added for debug
        print("→ JSON chain_names:", subunit_info['chain_names'])
        print("→ JSON sequence length:", len(subunit_info['sequence']))
        print("→ JSON sequence start:", subunit_info['sequence'][:20], "…")
        # Get node data (SubunitInfo object)
        iupred_shift = name_mapping[base_name]['start']-1
        start = node_data.start + iupred_shift
        end = node_data.end + iupred_shift
        sequence = node_data.sequence
        # Verify sequence matches the original
        #expected_sequence = subunit_info['sequence'][start-1:end]
        # Extract what we think is the “correct” slice
        full_seq = subunit_info['sequence']
        # If your SubunitInfo.start/end are *inclusive* 1-based, this is right:
        expected_sequence = full_seq[start - 1:end]
        # But if end were already exclusive, you’d want full_seq[start-1:end-1]…


        if not sequences_match(sequence, expected_sequence):
            print("==== DEBUG SEQUENCE MISMATCH ====")
            print(f"Node:           {node}")
            print(f"Start, end:     {start}, {end}")
            print(f"SubunitInfo.seq:   {sequence!r}  (len={len(sequence)})")
            print(f"Expected slice:    {expected_sequence!r}  (len={len(expected_sequence)})")
            print(f"Full sequence len: {len(full_seq)}")
            print(f"Chain names:       {subunit_info['chain_names']!r}")
            raise ValueError("Stopping here for debug")
        else:
            sequence = expected_sequence
        # if sequence != expected_sequence:
        #     raise ValueError(f"Sequence mismatch in node '{node}'. Expected: {expected_sequence}, Found: {sequence}")

        # Store high segment information
        node_index = len(high_segments_by_subunit[original_name]) + 1 #todo: changed here
        high_name = f"{subunit_info['name']}_high_{node_index}"
        high_segments_by_subunit[original_name].append({
            'node_name': node,
            'high_name': high_name,
            'start': start,
            'end': end
        })
        unified_subunits[high_name] = {
            'name': high_name,
            'sequence': sequence,
            'chain_names': subunit_info['chain_names'],
            'start_res': start
        }


    # Process each original subunit to create low segments
    for original_name, segments in high_segments_by_subunit.items():
        subunit_info = subunits_info[original_name]
        full_sequence = subunit_info['sequence']
        total_length = len(full_sequence)

        # Sort segments by start position
        segments.sort(key=lambda x: x['start'])

        # Track the end of the last processed segment
        last_end = 0
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
                low_name = f"{subunit_info['name']}_low_{low_segment_index}"

                # Extract sequence for low segment
                low_sequence = full_sequence[low_start-1:low_end]  # +1 because end is inclusive

                # Create low segment entry
                unified_subunits[low_name] = {
                    'name': low_name,
                    'sequence': low_sequence,
                    'chain_names': subunit_info['chain_names'],
                    'start_res': low_start  # Convert to 1-based indexing
                }

                low_segment_index += 1

            last_end = segment['end']
            print(f"Updated last_end to {last_end} after segment {segment}")

        # Check for gap after last high segment
        remaining_len = total_length - last_end

        if remaining_len > 0:
            # THRESHOLD CHECK: If the tail is tiny (e.g., < 5 residues), merge it backward
            if remaining_len < 5 and segments:
                # 1. Identify the last high segment we just created
                last_segment_dict = segments[-1]
                last_high_name = last_segment_dict['high_name']
                print(f"→ Merging tiny tail (size {remaining_len}) into {last_high_name}")

                # 2. Append the orphan sequence to it
                orphan_seq = full_sequence[last_end:]  # from last_end to end
                unified_subunits[last_high_name]['sequence'] += orphan_seq

            else:
                # Normal behavior: Create a new low-confidence segment
                low_start = last_end + 1
                low_end = total_length
                low_name = f"{subunit_info['name']}_low_{low_segment_index}"
                # Extract sequence for low segment
                low_sequence = full_sequence[low_start - 1:low_end]

                # Create low segment entry
                print(f"Creating low segment: {low_name} ({low_start}-{low_end})")
                unified_subunits[low_name] = {
                    'name': low_name,
                    'sequence': low_sequence,
                    'chain_names': subunit_info['chain_names'],
                    'start_res': low_start
                }

    sorted_unified_subunits = dict(
        sorted(
            unified_subunits.items(),
            key=lambda item: (item[1]['name'].split('_')[0], item[1]['start_res'])
        )
    )
    #output_folder = os.path.join(os.path.dirname(folder), 'combfold')
    os.makedirs(output_folder, exist_ok=True)
    output_json_path = os.path.join(output_folder, 'subunits_info.json')

    with open(output_json_path, 'w') as f:
        json.dump(sorted_unified_subunits, f, indent=4)


def save_connected_subunits_only(graph: nx.Graph, name_mapping: dict, subunits_info: dict,
                                 folder: str, output_subdir: str = 'combfold') -> None:
    """
    Create subunits_info.json containing ONLY subunits that have non-internal edges.
    Non-internal edge = edge between different chains (u[0] != v[0])

    Parameters
    ----------
    graph : nx.Graph
        Graph with SubunitInfo objects as node data
    name_mapping : dict
        Mapping of base names to chain names
    subunits_info : dict
        Dictionary containing subunit information
    folder : str
        Path to output folder
    output_subdir: str
        Subdirectory name for output
    """
    # Step 1: Identify nodes with NON-INTERNAL edges only
    connected_nodes = set()
    non_internal_edges = [
        (u, v) for u, v in graph.edges()
        if u[0] != v[0]  # Different chain prefixes
    ]

    print(f"\nTotal edges in graph: {len(graph.edges())}")
    print(f"Non-internal edges: {len(non_internal_edges)}")

    for edge in non_internal_edges:
        connected_nodes.add(edge[0])
        connected_nodes.add(edge[1])

    print(
        f"Found {len(connected_nodes)} connected nodes (via non-internal edges) out of {len(graph.nodes)} total nodes")
    print(f"Connected nodes: {sorted(connected_nodes)}")

    if not connected_nodes:
        print("WARNING: No non-internal edges found! Output file will be empty.")

    # Step 2: Build unified_subunits (same as save_subunits_info but filtered)
    unified_subunits = {}
    high_segments_by_subunit = defaultdict(list)

    # Only process connected nodes
    for node in connected_nodes:
        base_name = node.split('_')[0]
        node_data = graph.nodes[node]['data']
        original_name, subunit_info = find_original_subunit_info(
            base_name, name_mapping, subunits_info, node_data.start, node_data.end
        )

        iupred_shift = name_mapping[base_name]['start'] - 1
        start = node_data.start + iupred_shift
        end = node_data.end + iupred_shift
        sequence = node_data.sequence

        # Extract expected sequence
        full_seq = subunit_info['sequence']
        expected_sequence = full_seq[start - 1:end]

        if not sequences_match(sequence, expected_sequence):
            print(f"Warning: Sequence mismatch for node {node}")
            print(f"  Merged: {sequence[:50]}...")
            print(f"  Expected: {expected_sequence[:50]}...")
        else:
            sequence = expected_sequence

        # Store high segment
        node_index = len(high_segments_by_subunit[original_name]) + 1
        high_name = f"{subunit_info['name']}_high_{node_index}"
        high_segments_by_subunit[original_name].append({
            'node_name': node,
            'high_name': high_name,
            'start': start,
            'end': end
        })

        unified_subunits[high_name] = {
            'name': high_name,
            'sequence': sequence,
            'chain_names': subunit_info['chain_names'],
            'start_res': start
        }

    # # Step 3: Create low segments for gaps ONLY between connected high segments
    # for original_name, segments in high_segments_by_subunit.items():
    #     subunit_info = subunits_info[original_name]
    #     full_sequence = subunit_info['sequence']
    #     total_length = len(full_sequence)
    #
    #     segments.sort(key=lambda x: x['start'])
    #
    #     last_end = 0
    #     low_segment_index = 1
    #
    #     for segment in segments:
    #         current_start = segment['start']
    #
    #         # Gap before this segment
    #         if current_start > last_end + 1:
    #             low_start = last_end + 1
    #             low_end = current_start - 1
    #             low_name = f"{subunit_info['name']}_low_{low_segment_index}"
    #             low_sequence = full_sequence[low_start - 1:low_end]
    #
    #             unified_subunits[low_name] = {
    #                 'name': low_name,
    #                 'sequence': low_sequence,
    #                 'chain_names': subunit_info['chain_names'],
    #                 'start_res': low_start
    #             }
    #             low_segment_index += 1
    #
    #         last_end = segment['end']
    #
    #     # Gap after last segment
    #     remaining_len = total_length - last_end
    #     if remaining_len > 0:
    #         if remaining_len < 5 and segments:
    #             # Merge tiny tail into last high segment
    #             last_segment_dict = segments[-1]
    #             last_high_name = last_segment_dict['high_name']
    #             orphan_seq = full_sequence[last_end:]
    #             unified_subunits[last_high_name]['sequence'] += orphan_seq
    #         else:
    #             # Create new low segment
    #             low_start = last_end + 1
    #             low_end = total_length
    #             low_name = f"{subunit_info['name']}_low_{low_segment_index}"
    #             low_sequence = full_sequence[low_start - 1:low_end]
    #
    #             unified_subunits[low_name] = {
    #                 'name': low_name,
    #                 'sequence': low_sequence,
    #                 'chain_names': subunit_info['chain_names'],
    #                 'start_res': low_start
    #             }

    # Step 4: Sort and save
    sorted_unified_subunits = dict(
        sorted(
            unified_subunits.items(),
            key=lambda item: (item[1]['name'].split('_')[0], item[1]['start_res'])
        )
    )

    output_folder = os.path.join(os.path.dirname(folder), output_subdir)
    os.makedirs(output_folder, exist_ok=True)
    output_json_path = os.path.join(output_folder, 'subunits_info_connected_only.json')

    with open(output_json_path, 'w') as f:
        json.dump(sorted_unified_subunits, f, indent=4)

    print(f"\nSaved connected subunits (non-internal edges only) to: {output_json_path}")
    print(f"Total segments in filtered output: {len(sorted_unified_subunits)}")

    # Summary stats
    high_count = sum(1 for k in sorted_unified_subunits if '_high_' in k)
    low_count = sum(1 for k in sorted_unified_subunits if '_low_' in k)
    print(f"  High segments: {high_count}")
    print(f"  Low segments: {low_count}")

    # Show which chains are represented
    chains_represented = set()
    for info in sorted_unified_subunits.values():
        chains_represented.update(info['chain_names'])
    print(f"  Chains with interactions: {sorted(chains_represented)}")


def rename_graph_nodes (graph_to_rename, name_mapping):
    mapping_dict = {name:name_mapping[name]['chain_id'] for name in name_mapping}
    graph_to_rename = nx.relabel_nodes(graph_to_rename, mapping_dict)
    return graph_to_rename

if __name__ == "__main__":
    if len(sys.argv) >= 2:
        folder_path = os.path.abspath(sys.argv[1])

        # Check for optional "ipsae" flag
        use_ipsae_mode = False
        if "ipsae" in sys.argv:
            use_ipsae_mode = True
            # Remove "ipsae" from argv so it doesn't mess up file path indices if you used them loosely
            sys.argv.remove("ipsae")

        # Define Output Subdirectory based on mode
        if use_ipsae_mode:
            output_subdir_name = "combfold_ipsae"
            print(">>> RUNNING IN ipSAE MODE <<<")
        else:
            output_subdir_name = "combfold"
            print(">>> RUNNING IN pLDDT (Standard) MODE <<<")

        # (Standard loading logic...)
        mapping_path = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.join(os.path.dirname(folder_path),
                                                                                           'chain_id_mapping.json')
        original_subunits_path = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.path.join(
            os.path.dirname(folder_path), 'subunits_info.json')

        try:
            with open(mapping_path, 'r') as f:
                name_mapping = json.load(f)
            with open(original_subunits_path, 'r') as f:
                subunits_info = json.load(f)
        except Exception as e:
            print(f"Error loading JSONs: {e}")
            sys.exit(1)

        # 1. Create graphs (Passing the mode flag)
        graphs = create_graphs_from_folder(folder_path, use_ipsae=use_ipsae_mode)

        if not graphs:
            print("No graphs created. Exiting.")
            sys.exit(0)

        # 2. Merge graphs
        merged_graph = merge_graphs(graphs, name_mapping, subunits_info)

        # 3. Save with custom folder name
        save_subunits_info(merged_graph, name_mapping, subunits_info, folder_path, output_subdir=output_subdir_name)
        save_connected_subunits_only(merged_graph, name_mapping, subunits_info, folder_path, output_subdir_name)

        # 4. Check results (Look in the new folder)
        check_subunit_sequence_reconstruction(
            original_subunits_path,
            os.path.join(os.path.dirname(folder_path), output_subdir_name, "subunits_info.json")
        )

        final_graph = rename_graph_nodes(merged_graph, name_mapping)

        output_path = os.path.join(os.path.dirname(folder_path), output_subdir_name)
        pickle_file = os.path.join(output_path, "graph.pkl")
        pickle_file2 = os.path.join(output_path, "graph_before_renaming.pkl")

        with open(pickle_file, "wb") as f:
            pickle.dump(final_graph, f)

        with open(pickle_file2, "wb") as f:
            pickle.dump(merged_graph, f)

        show_circle(final_graph, os.path.dirname(folder_path))


        new_subunits_info_path = os.path.join(os.path.dirname(folder_path), output_subdir_name, "subunits_info.json")
        with open(new_subunits_info_path, "r") as f:
            new_subunits_info = json.load(f)
        show_subunit_circle_graph(final_graph, new_subunits_info, os.path.join(os.path.dirname(folder_path), output_subdir_name), mapping_path)

    else:
        print("usage: python create_graph_from_complex.py <folder_path> [mapping_json] [subunits_json] [ipsae]")