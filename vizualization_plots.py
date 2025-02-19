import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
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