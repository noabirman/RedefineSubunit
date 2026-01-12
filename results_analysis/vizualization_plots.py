import json
import os
import sys

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from matplotlib.patches import Wedge, Arc
from matplotlib import cm

from collections import defaultdict
import matplotlib.colors as mcolors
import math
import pickle

def spaced_colors(n):
    phi = 0.61803398875  # golden ratio conjugate
    return [(i * phi) % 1.0 for i in range(n)]


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

# def show_circle(graph: nx.Graph, folder_path: str):
#     """
#     Displays a graph in a circular layout with straight edges.
#
#     Args:
#         graph (nx.Graph): The graph to display.
#         folder_path (str): Path to save the output image.
#     """
#     plt.figure(figsize=(12, 12))
#
#     # Sort nodes alphabetically
#     sorted_nodes = sorted(graph.nodes())
#
#     # Generate positions for the nodes in circular layout
#     pos = nx.circular_layout(sorted_nodes)
#
#     # Draw nodes and labels
#     nx.draw_networkx_nodes(graph, pos, node_color='lightblue', node_size=300)
#     nx.draw_networkx_labels(graph, pos, font_size=5)
#
#     # Draw straight edges
#     nx.draw_networkx_edges(graph, pos, edge_color='gray')
#
#     # Hide axes
#     plt.axis('off')
#
#     # Save the plot
#     plt.savefig(os.path.join(folder_path, "curved_graph.png"))
def show_circle(graph: nx.Graph, folder_path: str):
    """
    Displays a graph in a circular layout.
    Nodes are colored by the first letter of their name.
    """
    n = graph.number_of_nodes()

    plt.figure(figsize=(16, 16))

    # Sort nodes alphabetically
    nodes = sorted(graph.nodes())

    # Circular layout
    pos = nx.circular_layout(nodes)

    # --- Adaptive sizing ---
    node_size = max(50, 3000 / max(n, 1))
    font_size = max(3, 12 / math.log2(n + 2))

    # --- Group nodes by first letter ---
    groups = sorted({node[0] for node in nodes})

    # Color map (supports many groups)
    #cmap = cm.get_cmap("tab20", len(groups))
    cmap = cm.get_cmap("hsv")  # or "turbo"
    # Map each group to a color
    group_to_color = {
        group: cmap(i / len(groups))
        for i, group in enumerate(groups)
    }

    # Assign color per node
    node_colors = [group_to_color[node[0]] for node in nodes]

    # Draw nodes
    nx.draw_networkx_nodes(
        graph,
        pos,
        node_color=node_colors,
        node_size=node_size,
        linewidths=0.5,
        edgecolors="black",
        alpha=0.9
    )

    # Draw labels
    nx.draw_networkx_labels(
        graph,
        pos,
        font_size=font_size
    )

    # Draw edges
    nx.draw_networkx_edges(
        graph,
        pos,
        edge_color="gray",
        width=0.5,
        alpha=0.6
    )

    plt.axis("off")

    plt.savefig(
        os.path.join(folder_path, "curved_graph.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

def draw_curve(ax, angle1, angle2, r_inner, color="gray", alpha=0.5, lw=0.5, curvature=0.6):
    x0, y0 = r_inner * np.cos(angle1), r_inner * np.sin(angle1)
    x2, y2 = r_inner * np.cos(angle2), r_inner * np.sin(angle2)

    # control point
    cx, cy = curvature * 0 + (1 - curvature) * (x0 + x2)/2, curvature * 0 + (1 - curvature) * (y0 + y2)/2

    t = np.linspace(0, 1, 50)
    x = (1-t)**2 * x0 + 2*(1-t)*t*cx + t**2 * x2
    y = (1-t)**2 * y0 + 2*(1-t)*t*cy + t**2 * y2
    ax.plot(x, y, color=color, alpha=alpha, lw=lw)


def show_subunit_circle_graph(graph: nx.Graph, subunit_info: dict, output_path: str,path_to_chain_mapping_json: str):
    """
    Circular arc graph:
    - Each protein is an arc (size ∝ total length)
    - Subunits split arcs
    - Edges connect subunit positions
    """

    # ---------- Parse subunits ----------
    proteins = {}
    chain_to_subunit = {}

    with open(path_to_chain_mapping_json, 'r') as f:
        chain_mapping = json.load(f)

    # Create reverse mapping: actual_chain_id -> list of graph_chain_ids
    actual_to_graph_chains = {}
    for graph_chain, info in chain_mapping.items():
        actual_chain = info['chain_id']
        actual_to_graph_chains.setdefault(actual_chain, []).append(graph_chain)

    for sub_name, info in subunit_info.items():
        protein = sub_name.split("_")[0]
        length = len(info["sequence"])

        for actual_chain in info["chain_names"]:
            # Get all graph chain IDs that map to this actual chain
            graph_chains = actual_to_graph_chains.get(actual_chain, [actual_chain])

            for graph_chain in graph_chains:
                # Create mapping using GRAPH chain IDs
                chain_key = f"{graph_chain}_{sub_name.split('_', 1)[1]}"  # e.g., "R_high_1"
                chain_to_subunit[chain_key] = sub_name  # Maps to "O43561_high_1"

        proteins.setdefault(protein, []).append({
            "name": sub_name,
            "length": length,
            "start_res": info["start_res"]
        })

    # Total lengths
    protein_lengths = {
        p: sum(s["length"] for s in subs)
        for p, subs in proteins.items()
    }
    total_length = sum(protein_lengths.values())

    # ---------- Color per protein ----------
    if len(proteins) <= 10:
        cmap = cm.get_cmap("tab10")
        colors = [cmap(i) for i in range(len(proteins))]
    elif len(proteins) <= 20:
        cmap = cm.get_cmap("tab20")
        colors = [cmap(i) for i in range(len(proteins))]
    else:
        # For many proteins, use improved HSV
        golden_ratio = 0.618033988749895
        colors = [cm.get_cmap("hsv")((i * golden_ratio) % 1.0) for i in range(len(proteins))]

    protein_colors = {
        protein: colors[i]
        for i, protein in enumerate(sorted(proteins))
    }

    # ---------- Figure ----------
    fig, ax = plt.subplots(figsize=(28, 28))
    ax.set_aspect("equal")
    ax.axis("off")

    radius = 25.0
    width = 2.0

    angle_start = 0
    subunit_angles = {}  # subunit -> angle (for edges)

    # Gap between proteins (in degrees)
    gap_size = 2.0  # Adjust this value for larger/smaller gaps
    total_gap = gap_size * len(proteins)
    # Adjust total angle to account for gaps
    total_angle = 360 - total_gap

    # ---------- Draw arcs ----------
    for protein in sorted(proteins):
        p_len = protein_lengths[protein]
        p_angle = total_angle * p_len / total_length
        p_start = angle_start
        p_end = angle_start + p_angle


        # ---------- Subunit divisions ----------
        sub_start = p_start
        for sub in proteins[protein]:
            sub_angle = p_angle * sub["length"] / p_len

            # Set alpha for low/high
            if "low" in sub["name"]:
                alpha = 0.4
            else:
                alpha = 1.0

            # Draw the subunit wedge
            wedge = Wedge(
                center=(0, 0),
                r=radius,
                theta1=sub_start,
                theta2=sub_start + sub_angle,
                width=width,
                facecolor=protein_colors[protein],
                alpha=alpha,
                edgecolor="black",
                linewidth=1.5
            )
            ax.add_patch(wedge)

            sub_mid = sub_start + sub_angle / 2
            subunit_angles[sub["name"]] = math.radians(sub_mid)

            # Stronger divider line at START of each subunit
            theta = math.radians(sub_start)
            ax.plot(
                [radius * math.cos(theta), (radius - width) * math.cos(theta)],
                [radius * math.sin(theta), (radius - width) * math.sin(theta)],
                color="black",
                linewidth=2.0,
                alpha=1.0
            )
            sub_start += sub_angle

        # Add final divider at end of protein
        theta = math.radians(sub_start)
        ax.plot(
            [radius * math.cos(theta), (radius - width) * math.cos(theta)],
            [radius * math.sin(theta), (radius - width) * math.sin(theta)],
            color="black",
            linewidth=2.0,
            alpha=1.0
        )
        # ---------- Add residue ruler (OUTSIDE the arc) ----------
        # Draw ruler circle arc for this protein
        ruler_radius = radius + 0.5
        ruler_arc = Arc(
            xy=(0, 0),
            width=2 * ruler_radius,
            height=2 * ruler_radius,
            angle=0,
            theta1=p_start,
            theta2=p_end,
            color='black',
            linewidth=1.0,
            alpha=0.6
        )
        ax.add_patch(ruler_arc)

        # Determine tick interval based on protein length
        if p_len < 200:
            tick_interval = 50
        elif p_len < 500:
            tick_interval = 100
        elif p_len < 1000:
            tick_interval = 200
        else:
            tick_interval = 250

        # Draw ticks and labels OUTSIDE the arc
        for residue in range(0, p_len + 1, tick_interval):
            if residue > p_len:
                continue

            # Calculate angle for this residue
            res_angle = p_start + p_angle * (residue / p_len)
            theta = math.radians(res_angle)

            # Draw tick mark (pointing outward)
            tick_start_r = ruler_radius
            tick_end_r = ruler_radius + 0.8
            ax.plot(
                [tick_start_r * math.cos(theta), tick_end_r * math.cos(theta)],
                [tick_start_r * math.sin(theta), tick_end_r * math.sin(theta)],
                color="black",
                linewidth=1.0,
                alpha=0.8
            )

            # Add residue number label (outside the tick, parallel to ruler)
            label_r = ruler_radius + 1.5
            rotation = res_angle  # Parallel to the circle (tangent)

            # Adjust text rotation for readability (flip if upside down)
            if 90 < res_angle < 270:
                rotation = res_angle + 180  # Flip to keep text right-side up

            ax.text(
                label_r * math.cos(theta),
                label_r * math.sin(theta),
                str(residue),
                ha='center',
                va='center',
                fontsize=14,
                rotation=rotation,
                rotation_mode='anchor',
                color='black'
            )
        # Protein label
        mid_angle = math.radians((p_start + p_end) / 2)
        label_radius = ruler_radius + 4.0  # Further out, past the ruler
        ax.text(
            label_radius * math.cos(mid_angle),
            label_radius * math.sin(mid_angle),
            protein,
            ha="center",
            va="center",
            fontsize=18,
            fontweight='bold'
        )
        angle_start += p_angle + gap_size

    # ---------- Draw edges ----------
    # ---------- Better Edge Debug ----------
    print("\n=== Edge Debug Info ===")
    print(f"Total edges in graph: {len(graph.edges())}")
    print(f"Subunit angles available: {len(subunit_angles)}")
    print(f"chain_to_subunit mappings: {len(chain_to_subunit)}")

    # Sample data
    print(f"\nSample graph edges (first 5): {list(graph.edges())[:5]}")
    print(f"Sample chain_to_subunit keys (first 10): {list(chain_to_subunit.keys())[:10]}")
    print(f"Sample subunit_angles keys (first 10): {list(subunit_angles.keys())[:10]}")


    graph_chains = set()
    count_edges = 0
    for u, v in graph.edges():
        su = chain_to_subunit.get(u)
        sv = chain_to_subunit.get(v)

        graph_chains.add(u)
        graph_chains.add(v)


        if su not in subunit_angles or sv not in subunit_angles:
            continue

        a1 = subunit_angles[su]
        a2 = subunit_angles[sv]

        # Determine if edge is inter-protein or intra-protein
        protein_u = su.split("_")[0]  # Extract protein from subunit name
        protein_v = sv.split("_")[0]

        if protein_u == protein_v:
            # Intra-protein edge (within same protein)
            edge_color = "gray"
            edge_alpha = 0.6
            #edge_width = 0.5
        else:
            # Inter-protein edge (between different proteins)
            edge_color = "magenta"  # or "purple", "blue", etc.
            edge_alpha = 0.9
            count_edges+=1
            #edge_width = 1.0

        draw_curve(ax, a1, a2, r_inner=radius - width / 2,
                   color=edge_color, alpha=edge_alpha, lw=0.5, curvature=0.3)


    # Find which ones are missing from chain_to_subunit
    missing_chains = graph_chains - set(chain_to_subunit.keys())

    print(f"\nTotal unique chains in graph: {len(graph_chains)}")
    print(f"Chains in chain_to_subunit: {len(chain_to_subunit)}")
    print(f"Missing chains ({len(missing_chains)}): {sorted(missing_chains)[:20]}")
    print(f"number of non internal edges: {count_edges}")
    print("======================\n")
    # ---------- Save ----------
    plt.savefig(os.path.join(output_path, "subunit_graph.png"), dpi=450, bbox_inches="tight")  # Even higher DPI
    plt.close()


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


def group_nodes_by_chain_name(G: nx.Graph, nodes_data):
    # Group nodes by chain name
    chain_groups = defaultdict(list)
    for node, attrs in nodes_data.items():
        chain = attrs['data'].chain_names[0]  # assuming one chain per subunit
        start = attrs['data'].start
        chain_groups[chain].append((node, start))

    # Sort each group by start position
    for chain in chain_groups:
        chain_groups[chain] = sorted(chain_groups[chain], key=lambda x: x[1])
    return chain_groups


def plot_graph_by_chain(G: nx.Graph):
    #chain_groups = group_nodes_by_chain_name(G, nodes_data)

    # === Group nodes by chain and sort by start position ===
    chain_groups = defaultdict(list)
    for node in G.nodes(data=True):
        chain = node[1]['data'].chain_names[0]
        start = node[1]['data'].start
        chain_groups[chain].append((node[0], start))

    # Sort nodes within each chain by start position
    for chain in chain_groups:
        chain_groups[chain] = sorted(chain_groups[chain], key=lambda x: x[1])

    # Assign positions to each node
    pos = {}
    x_spacing = 5
    y_spacing = 5
    x = 0
    for i, (chain, nodes) in enumerate(chain_groups.items()):
        y = 0
        for node_id, _ in nodes:
            pos[node_id] = (x, y)
            y -= y_spacing
        x += x_spacing

    # === Assign colors to each chain ===
    chains = sorted(chain_groups.keys())
    color_map = cm.get_cmap('tab10', len(chains))
    chain_to_color = {chain: mcolors.to_hex(color_map(i)) for i, chain in enumerate(chains)}

    node_colors = []
    for node in G.nodes():
        chain = G.nodes[node]['data'].chain_names[0]
        node_colors.append(chain_to_color[chain])

    # === Draw the graph ===
    plt.figure(figsize=(12, 8), constrained_layout=True)
    nx.draw(
        G, pos, with_labels=True, node_size=600,
        font_size=8, node_color=node_colors, edgecolors='black'
    )

    # === Legend ===
    legend_handles = [Patch(color=color, label=chain) for chain, color in chain_to_color.items()]
    plt.legend(handles=legend_handles, title='Chains', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.title("Subunits arranged and colored by chain")
    plt.axis("off")
    #plt.tight_layout()
    plt.savefig("G_by_chain.png")

if __name__ == "__main__":
    # test the function
    folder = os.path.abspath(sys.argv[1])
    parent_dir = os.path.dirname(folder)
    chain_mapping_path = os.path.join(parent_dir, 'chain_id_mapping.json')

    subunits_info_path = os.path.join(folder, "subunits_info.json")

    with open(subunits_info_path, 'r') as f:
        subunits_info = json.load(f)


    with open(os.path.join(folder,"graph.pkl"), "rb") as f:
        graph = pickle.load(f)


    show_subunit_circle_graph(graph, subunits_info, folder,chain_mapping_path)

