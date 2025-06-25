#!/usr/bin/env python3
"""
Simple Plotting Script for OpenSees Cantilever Beam Results
This script creates visualizations using matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.patches as patches

def parse_nodes_file(filename):
    """Parse the nodes.txt file to extract coordinates"""
    nodes = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('Node:'):
            # Extract node ID
            parts = line.split()
            node_id = int(parts[1])
            
            # Look for coordinates in next line
            if i + 1 < len(lines):
                coord_line = lines[i + 1].strip()
                if 'Coordinates' in coord_line:
                    # Extract coordinates
                    coord_start = coord_line.find(':') + 1
                    coord_str = coord_line[coord_start:].strip()
                    coords = [float(x) for x in coord_str.split()]
                    if len(coords) >= 2:
                        nodes.append((node_id, coords[0], coords[1]))
        i += 1
    
    return nodes

def parse_elements_file(filename):
    """Parse the elements.txt file to extract element connectivity"""
    elements = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('FourNodeQuad, element id:'):
            # Extract element ID
            parts = line.split()
            ele_id = int(parts[3])
            
            # Look for connected nodes in next line
            if i + 1 < len(lines):
                nodes_line = lines[i + 1].strip()
                if 'Connected external nodes:' in nodes_line:
                    # Extract node IDs
                    nodes_start = nodes_line.find(':') + 1
                    nodes_str = nodes_line[nodes_start:].strip()
                    node_ids = [int(x) for x in nodes_str.split()]
                    elements.append((ele_id, node_ids))
        i += 1
    
    return elements

def parse_displacement_file(filename):
    """Parse the displacements.out file to extract displacement data"""
    with open(filename, 'r') as f:
        line = f.readline().strip()
    
    # Split the line into values
    values = [float(x) for x in line.split()]
    
    # First value is time step, rest are displacements
    time_step = values[0]
    displacements = values[1:]
    
    return time_step, displacements

def parse_element_stress(filename):
    """Parse elements.txt to extract and average stress from Gauss points."""
    element_stresses = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        if 'Stress (xx yy xy)' in lines[i]:
            gauss_stresses = []
            # Next 4 lines should be Gauss point stresses
            for j in range(1, 5):
                if i + j < len(lines) and 'Gauss point' in lines[i+j]:
                    try:
                        parts = lines[i+j].split(':')[1].strip().split()
                        stress_vals = [float(p) for p in parts]
                        gauss_stresses.append(stress_vals)
                    except (ValueError, IndexError):
                        continue
            
            if len(gauss_stresses) > 0:
                # Average the stresses (xx, yy, xy)
                avg_stress = np.mean(gauss_stresses, axis=0)
                element_stresses.append(avg_stress)
            
            # Move index past this stress block
            i += 4
        else:
            i += 1
            
    return np.array(element_stresses)

def plot_displacement_history():
    """Plot displacement history from OpenSees output"""
    
    if not os.path.exists('RESULTS/displacements.out'):
        print("ERROR: displacements.out not found")
        return
    
    # Parse displacement data
    time_step, disp_data = parse_displacement_file('RESULTS/displacements.out')
    
    # Parse nodes to get coordinates
    nodes = parse_nodes_file('RESULTS/nodes.txt')
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Plot displacement for each node (every 2 values = x and y displacement)
    for i in range(0, min(20, len(disp_data)), 2):  # Plot first 10 nodes
        if i + 1 < len(disp_data):
            node_id = i // 2 + 1
            disp_x = disp_data[i]
            disp_y = disp_data[i + 1]
            plt.plot([time_step], [disp_y], 'o', label=f'Node {node_id} (Y)', markersize=8)
    
    plt.xlabel('Time Step')
    plt.ylabel('Displacement (inches)')
    plt.title('Displacement History')
    plt.legend()
    plt.grid(True)
    plt.savefig('RESULTS/displacement_history.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_mesh_and_deformed():
    """Plot original and deformed mesh with filled Q4 elements and connected nodes"""
    
    if not os.path.exists('RESULTS/nodes.txt'):
        print("ERROR: nodes.txt not found")
        return
    
    if not os.path.exists('RESULTS/elements.txt'):
        print("ERROR: elements.txt not found")
        return
    
    if not os.path.exists('RESULTS/displacements.out'):
        print("ERROR: displacements.out not found")
        return
    
    # Parse nodes, elements, and displacements
    nodes = parse_nodes_file('RESULTS/nodes.txt')
    elements = parse_elements_file('RESULTS/elements.txt')
    time_step, disp_data = parse_displacement_file('RESULTS/displacements.out')
    
    if len(nodes) == 0:
        print("ERROR: No nodes found in nodes.txt")
        return
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Create node dictionary for easy lookup
    node_dict = {node[0]: (node[1], node[2]) for node in nodes}
    
    # Plot original mesh
    # Plot nodes
    x_coords = [node[1] for node in nodes]
    y_coords = [node[2] for node in nodes]
    ax1.scatter(x_coords, y_coords, c='black', s=50, zorder=5)
    
    # Plot elements as filled polygons
    for ele_id, node_ids in elements:
        if len(node_ids) == 4:  # Q4 element
            # Get coordinates of the 4 nodes
            coords = []
            for node_id in node_ids:
                if node_id in node_dict:
                    coords.append(node_dict[node_id])
            
            if len(coords) == 4:
                # Create polygon
                poly_x = [coord[0] for coord in coords]
                poly_y = [coord[1] for coord in coords]
                
                # Fill the element
                ax1.fill(poly_x, poly_y, alpha=0.3, color='lightblue', edgecolor='blue', linewidth=1)
                
                # Connect nodes with lines
                for i in range(4):
                    x1, y1 = coords[i]
                    x2, y2 = coords[(i + 1) % 4]
                    ax1.plot([x1, x2], [y1, y2], 'b-', linewidth=1, zorder=3)
    
    ax1.set_title('Original Mesh')
    ax1.set_xlabel('X (inches)')
    ax1.set_ylabel('Y (inches)')
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')
    
    # Plot deformed mesh
    # Calculate deformed coordinates
    deformed_coords = {}
    for i, node in enumerate(nodes):
        if 2*i < len(disp_data):
            def_x = node[1] + disp_data[2*i]
            def_y = node[2] + disp_data[2*i + 1]
        else:
            def_x = node[1]
            def_y = node[2]
        deformed_coords[node[0]] = (def_x, def_y)
    
    # Plot deformed nodes
    def_x_coords = [deformed_coords[node[0]][0] for node in nodes]
    def_y_coords = [deformed_coords[node[0]][1] for node in nodes]
    ax2.scatter(def_x_coords, def_y_coords, c='black', s=50, zorder=5)
    
    # Plot deformed elements
    for ele_id, node_ids in elements:
        if len(node_ids) == 4:  # Q4 element
            # Get deformed coordinates of the 4 nodes
            coords = []
            for node_id in node_ids:
                if node_id in deformed_coords:
                    coords.append(deformed_coords[node_id])
            
            if len(coords) == 4:
                # Create polygon
                poly_x = [coord[0] for coord in coords]
                poly_y = [coord[1] for coord in coords]
                
                # Fill the element
                ax2.fill(poly_x, poly_y, alpha=0.3, color='lightcoral', edgecolor='red', linewidth=1)
                
                # Connect nodes with lines
                for i in range(4):
                    x1, y1 = coords[i]
                    x2, y2 = coords[(i + 1) % 4]
                    ax2.plot([x1, x2], [y1, y2], 'r-', linewidth=1, zorder=3)
    
    # Find tip node and annotate displacement magnitude
    tip_node = max(nodes, key=lambda x: x[1])  # Node with maximum x coordinate
    tip_node_id = tip_node[0]
    
    if tip_node_id in deformed_coords:
        tip_orig = node_dict[tip_node_id]
        tip_def = deformed_coords[tip_node_id]
        
        # Calculate displacement magnitude
        disp_mag = np.sqrt((tip_def[0] - tip_orig[0])**2 + (tip_def[1] - tip_orig[1])**2)
        
        # Annotate tip displacement
        ax2.annotate(f'Tip Displacement: {disp_mag:.4f} in', 
                    xy=tip_def, xytext=(tip_def[0] - 20, tip_def[1] + 5),
                    arrowprops=dict(arrowstyle='->', color='red', lw=2),
                    fontsize=10, color='red', weight='bold',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        
        # Mark tip node
        ax2.scatter([tip_def[0]], [tip_def[1]], c='red', s=100, zorder=6, marker='*')
    
    ax2.set_title('Deformed Mesh')
    ax2.set_xlabel('X (inches)')
    ax2.set_ylabel('Y (inches)')
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig('RESULTS/mesh_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_stress_distribution():
    """Plot sigma_xx, sigma_yy, and tau_xy stress distributions on undeformed mesh as subplots"""
    if not os.path.exists('RESULTS/elements.txt'):
        print("ERROR: elements.txt not found")
        return

    # Parse nodes, elements, and displacements
    nodes = parse_nodes_file('RESULTS/nodes.txt')
    elements = parse_elements_file('RESULTS/elements.txt')
    
    # Parse element stress
    stress = parse_element_stress('RESULTS/elements.txt')
    
    if stress.shape[0] == 0:
        print("ERROR: No stress data found in elements.txt. The analysis may have failed.")
        return
        
    if stress.shape[0] != len(elements):
        print(f"Warning: Number of stresses ({stress.shape[0]}) does not match number of elements ({len(elements)}).")
        return

    # Create node dictionary for easy lookup
    node_dict = {node[0]: (node[1], node[2]) for node in nodes}

    # Get axis limits from node coordinates
    all_x = [node[1] for node in nodes]
    all_y = [node[2] for node in nodes]
    xlim = (min(all_x), max(all_x))
    ylim = (min(all_y), max(all_y))

    # Prepare stress components and titles
    stress_names = [r'$\sigma_{xx}$', r'$\sigma_{yy}$', r'$\tau_{xy}$']
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True, sharey=True)
    for k, ax in enumerate(axes):
        cmap = plt.get_cmap('jet')
        vals = stress[:, k]
        vmin = np.min(vals)
        vmax = np.max(vals)
        max_idx = np.argmax(vals)
        min_idx = np.argmin(vals)
        centroids = []
        for idx, (ele_id, node_ids) in enumerate(elements):
            if len(node_ids) == 4 and idx < len(vals):
                coords = [node_dict[nid] for nid in node_ids if nid in node_dict]
                if len(coords) == 4:
                    color_val = (vals[idx] - vmin) / (vmax - vmin) if vmax > vmin else 0.5
                    poly_x = [c[0] for c in coords]
                    poly_y = [c[1] for c in coords]
                    ax.fill(poly_x, poly_y, color=cmap(color_val), alpha=0.8, edgecolor='k', linewidth=0.5)
                    # Store centroid for annotation
                    centroid = (np.mean(poly_x), np.mean(poly_y))
                    centroids.append(centroid)
                else:
                    centroids.append((None, None))
            else:
                centroids.append((None, None))
        ax.set_title(f'{stress_names[k]} Stress')
        ax.set_xlabel('X (inches)')
        if k == 0:
            ax.set_ylabel('Y (inches)')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        # Add colorbar with height matching the axis
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.1)
        cb = fig.colorbar(sm, cax=cax, label=stress_names[k] + ' (ksi)')
        # Annotate max and min
        if centroids[max_idx][0] is not None:
            ax.annotate(f"max: {vmax:.2f}", xy=centroids[max_idx], xytext=(centroids[max_idx][0], centroids[max_idx][1]),
                        color='red', fontsize=10, fontweight='bold', ha='center', va='center',
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="red", lw=1, alpha=0.7))
        if centroids[min_idx][0] is not None:
            ax.annotate(f"min: {vmin:.2f}", xy=centroids[min_idx], xytext=(centroids[min_idx][0], centroids[min_idx][1]),
                        color='blue', fontsize=10, fontweight='bold', ha='center', va='center',
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="blue", lw=1, alpha=0.7))
    plt.tight_layout()
    plt.savefig('RESULTS/stress_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    """Main function"""
    print("=== OpenSees Results Visualization ===")
    
    # Check if RESULTS directory exists
    if not os.path.exists('RESULTS'):
        print("ERROR: RESULTS directory not found. Run the OpenSees analysis first.")
        return
    
    # Create visualizations
    print("Creating displacement history plot...")
    plot_displacement_history()
    
    print("Creating mesh comparison plot...")
    plot_mesh_and_deformed()
    
    print("Creating stress distribution plot...")
    plot_stress_distribution()
    
    print("\nVisualization complete!")
    print("Check the RESULTS directory for output images:")
    print("  - displacement_history.png")
    print("  - mesh_comparison.png")
    print("  - stress_distribution.png")

if __name__ == "__main__":
    main() 