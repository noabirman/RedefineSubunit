import os

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from collections import defaultdict
import matplotlib.cm as cm
import matplotlib.colors as mcolors

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