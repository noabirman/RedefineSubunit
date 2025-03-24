import os
import sys
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
from complex_graph import graph, SubunitName, SubunitInfo
from typing import List


def overlap(v1: SubunitInfo, v2: SubunitInfo) -> bool:
    """Check if two SubunitInfo nodes overlap in at least one chain."""
    return any(chain in v2.chain_names for chain in v1.chain_names) and not (v1.end < v2.start or v2.end < v1.start)


def merge_graphs(graphs: List[nx.Graph]) -> nx.Graph:
    """Merge nodes with overlapping indices in the same chain across multiple graphs."""
    merged_graph = nx.Graph()

    # Step 1: Build an overlap graph
    overlap_graph = nx.Graph()  # Temporary graph for finding connected components

    # Add all nodes to the overlap graph
    all_nodes = {}
    for G in graphs:
        for node, data in G.nodes(data=True):
            all_nodes[node] = data['data']
            overlap_graph.add_node(node, data=data['data'])

    # Add edges between overlapping nodes
    nodes_list = list(all_nodes.items())  # Convert to list for pairwise comparison
    for i in range(len(nodes_list)):
        name1, subunit1 = nodes_list[i]
        for j in range(i + 1, len(nodes_list)):
            name2, subunit2 = nodes_list[j]
            if overlap(subunit1, subunit2):
                overlap_graph.add_edge(name1, name2)

    # Step 2: Find connected components (groups of merged nodes)
    connected_components = list(nx.connected_components(overlap_graph))

    # Step 3: Merge nodes in each component
    merged_nodes = {}
    node_mapping = {}

    for component in connected_components:
        subunits = [all_nodes[name] for name in component]

        # Merge properties
        merged_name = SubunitName(f"{subunits[0].name[0:-1]}: {subunits[0].start}-{subunits[0].end}")
        merged_chains = sorted(set(chain for subunit in subunits for chain in subunit.chain_names))
        merged_start = min(subunit.start for subunit in subunits)
        merged_end = max(subunit.end for subunit in subunits)
        # merged_sequence = "".join(subunit.sequence for subunit in subunits)  # Optional

        merged_subunit = SubunitInfo(
            name=merged_name,
            chain_names=merged_chains,
            start=merged_start,
            end=merged_end,
            sequence=''
        )

        merged_nodes[merged_name] = merged_subunit

        # Map original nodes to merged node name
        for name in component:
            node_mapping[name] = merged_name

    # Step 4: Add merged nodes to the final graph
    for merged_name, subunit in merged_nodes.items():
        merged_graph.add_node(merged_name, data=subunit)

    # Step 5: Update edges based on merged nodes
    for G in graphs:
        for u, v in G.edges():
            merged_u = node_mapping[u]
            merged_v = node_mapping[v]
            merged_graph.add_edge(merged_u, merged_v)

    return merged_graph

def create_graphs_from_folder(folder):
    graphs = []
    # Iterate over all folders in the folder
    for root, dirs, files in os.walk(folder):
        for dir in dirs:
            data_path = os.path.join(root, dir, f"{dir}_confidences.json")
            data_path = os.path.abspath(data_path)
            structure_path = os.path.join(root, dir, f"{dir}_model.cif")
            structure_path = os.path.abspath(structure_path)
            graph1 = graph(structure_path, data_path, af_version = '3')
            graphs.append(graph1)
    return graphs


def show_graph(graph: nx.Graph):
    """Print the nodes and edges of a graph."""
    print("Nodes:")
    for node, data in graph.nodes(data=True):
        print(f"  {node}: {data['data']}")

    print("Edges:")
    for u, v in graph.edges():
        print(f"  {u} -- {v}")
    # Draw the graph
    plt.figure(figsize=(6, 6))
    nx.draw(graph, with_labels=True, node_color='lightblue', edge_color='gray', node_size=800, font_size=10)
    plt.show()





# main
if __name__ == "__main__":
    if len(sys.argv) == 2:
        folder_name = os.path.abspath(sys.argv[1])
        graphs = create_graphs_from_folder(folder_name)
        merged_graph = merge_graphs(graphs)
        show_graph(merged_graph)
    else:
        print("usage: <script> enter folder_name")


