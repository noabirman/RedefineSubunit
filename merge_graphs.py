import os
import sys
import random
from collections import defaultdict
import networkx as nx
from networkx.readwrite import json_graph
import matplotlib.pyplot as plt
from complex_graph import graph, SubunitName, SubunitInfo
from typing import List
import numpy as np


def overlap(v1: SubunitInfo, v2: SubunitInfo) -> bool:
    """Check if two SubunitInfo nodes overlap in at least one chain."""
    return any(chain in v2.chain_names for chain in v1.chain_names) and not (v1.end < v2.start or v2.end < v1.start)


def merge_graphs(graphs: List[nx.Graph]) -> nx.Graph:
    print("Merging graphs")
    """Merge nodes with overlapping indices in the same chain across multiple graphs."""
    merged_graph = nx.Graph()

    # Step 1: Build an overlap graph
    print("Building overlap graph")
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
    print("Finding connected components")
    connected_components = list(nx.connected_components(overlap_graph))

    # Step 3: Merge nodes in each component
    print("Merging nodes")
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
    print(f"Creating graphs from folder {folder}")
    graphs = []
    # Iterate over all folders in the folder
    for item in os.listdir(folder):
        item_path = os.path.join(folder, item)
        if os.path.isdir(item_path):
            print(f"Processing {item}")
            data_path = os.path.join(item_path, f"{item}_confidences.json")
            data_path = os.path.abspath(data_path)
            structure_path = os.path.join(item_path, f"{item}_model.cif")
            structure_path = os.path.abspath(structure_path)
            graph1 = graph(structure_path, data_path, af_version = '3')
            graphs.append(graph1)
    return graphs

def show_graph_with_spacing(graph: nx.Graph, folder_path: str, name:str):
    """
    Displays the graph with more space between nodes using spring_layout and adjusting the k parameter.

    Args:
        graph (nx.Graph): The graph to display.
        folder_path (str): Path to save the output image.
    """
    # Calculate the figure size dynamically based on the number of nodes
    num_nodes = len(graph.nodes())
    figsize = (max(6, num_nodes / 10), max(6, num_nodes / 10))  # Adjust figure size

    # Dynamically calculate proportional node size and font size
    node_size = max(100, 500 / num_nodes)
    font_size = max(8, 20 / num_nodes)

    # Use spring_layout to position nodes with more space (increase k value)
    pos = nx.spring_layout(graph, k=0.3, iterations=50)  # Increase 'k' for more space between nodes

    # Draw the graph with adjusted layout and sizes
    plt.figure(figsize=figsize)
    nx.draw(graph, pos, with_labels=True, node_size=node_size, font_size=font_size, edge_color='gray',
            node_color='lightblue')

    # Save the graph image
    save_path = os.path.join(folder_path, name+"graph_with_spacing.png")
    plt.savefig(save_path, bbox_inches='tight')
    plt.show()

def show_circle(graph: nx.Graph, folder_path:str, name:str):
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch

    # Sample graph creation for demonstration
    graph = nx.Graph()
    graph.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'A'), ('A', 'D'), ('B', 'D'), ('C', 'D')])

    plt.figure(figsize=(12, 12))

    # Sort nodes alphabetically
    sorted_nodes = sorted(graph.nodes())

    # Generate positions for the nodes in circular layout
    pos = nx.circular_layout(sorted_nodes)

    # Draw nodes
    nx.draw_networkx_nodes(graph, pos, node_color='lightblue', node_size=300)
    nx.draw_networkx_labels(graph, pos, font_size=5)

    # Function to compute the control points for a cubic Bezier curve
    def compute_control_points(p1, p2, curve_factor):
        midpoint = (p1 + p2) / 2
        dir_v = p2 - p1
        perp_v = np.array([dir_v[1], -dir_v[0]])
        perp_v /= np.linalg.norm(perp_v)
        contrl1 = midpoint - curve_factor * perp_v
        contrl2 = midpoint + curve_factor * perp_v
        return contrl1, contrl2

    # Draw edges with dynamic curvature
    for edge in graph.edges():
        p1 = np.array(pos[edge[0]])
        p2 = np.array(pos[edge[1]])

        # Distance between nodes
        distance = np.linalg.norm(p1 - p2)

        # Curve factor: closer nodes = more curve
        max_dist = 2.0  # Max distance for unit circle nodes
        normalized_dist = distance / max_dist
        curve_factor = 0.5 * (1 - normalized_dist)

        # Compute control points for cubic Bezier curve
        contrl1, contrl2 = compute_control_points(p1, p2, curve_factor)

        # Create Bezier path
        verts = [p1, contrl1, contrl2, p2]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path = Path(verts, codes)

        # Add path to plot
        patch = PathPatch(path, edgecolor='gray', fill=False)
        plt.gca().add_patch(patch)

    # Hide axes
    plt.axis('off')
    plt.savefig(os.path.join(folder_path,f"curved_{name}.png"))

def show_graph(graph: nx.Graph, folder_path:str):
    """Print the nodes and edges of a graph."""
    print(f"total number of Nodes: {graph.number_of_nodes()}")
    print(f"total number of Edges: {graph.number_of_edges()}")

    print("Nodes:")
    for node, data in graph.nodes(data=True):
        print(f"  {node}: {data['data']}")

    print("Edges:")
    for u, v in graph.edges():
        print(f"  {u} -- {v}")

    # # Extract chain names from node attributes
    # chain_names = {graph.nodes[node][1]['chain'] for node in graph.nodes(data=True)}
    # # Assign unique colors to each chain
    # chain_colors = {chain: (random.random(), random.random(), random.random()) for chain in chain_names}
    # # Get node colors based on their chain
    # node_colors = [chain_colors[graph.nodes[node][1]['chain']] for node in graph.nodes(data=True)]

    # Draw the graph
    plt.figure(figsize=(12, 12))
    nx.draw(graph, with_labels=True, node_color='lightblue', edge_color='gray', node_size=300, font_size=5)
    plt.savefig(os.path.join(folder_path, "merged_graph.png"))
    plt.show()

    nodes_with_edges = [node for node, degree in graph.degree() if degree > 0]
    # Create a subgraph with only these nodes
    graph_filtered = graph.subgraph(nodes_with_edges).copy()
    plt.figure(figsize=(12, 12))
    nx.draw(graph_filtered, with_labels=True, node_color='lightblue', edge_color='gray', node_size=800, font_size=10)
    plt.savefig(os.path.join(folder_path, "edges_only_merged_graph.png"))
    plt.show() #gg
    show_graph_with_spacing(graph, folder_path,"reg_")
    show_graph_with_spacing(graph_filtered, folder_path, "filtered")

# main
if __name__ == "__main__":
    if len(sys.argv) == 2:
        folder_path = os.path.abspath(sys.argv[1])
        folder_name = os.path.split(folder_path)[1]
        graphs = create_graphs_from_folder(folder_path)
        merged_graph = merge_graphs(graphs)
        json_graph.dumps(merged_graph)
        nodes_with_edges = [node for node, degree in merged_graph.degree() if degree > 0]
        sub_graph = merged_graph.subgraph(nodes_with_edges).copy()
        show_graph(merged_graph, folder_path)
        show_circle(merged_graph, folder_path, 'graph')
        show_circle(sub_graph, folder_path, 'sub_graph')

    else:
        print("usage: <script> enter folder_name")


