#!/usr/bin/env python3
"""
OpenSees Earthquake Model Results Visualization
Plots displacement history, mesh comparison, stress distribution, and deformation animation
Updated for Q4 elements and corrected output format
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from matplotlib.patches import Polygon

def parse_nodes_file(filename):
    """Parse the eq_nodes.txt file to extract node coordinates"""
    if not os.path.exists(filename):
        print(f"ERROR: {filename} not found.")
        return []
    nodes = []
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith('Node'):
                parts = line.split()
                if len(parts) >= 4:
                    node_id = int(parts[1].replace(':', ''))
                    # Extract coordinates - new format: "Node 1: x y"
                    # Skip the first two parts ("Node" and "1:")
                    x = float(parts[2])
                    y = float(parts[3])
                    nodes.append((node_id, x, y))
    except Exception as e:
        print(f"Error reading nodes from {filename}: {e}")
    return nodes

def parse_elements_file(filename):
    """Parse the eq_elements.txt file to extract element connectivity"""
    if not os.path.exists(filename):
        print(f"ERROR: {filename} not found.")
        return []
    elements = []
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith('Element') and 'nodes' in line:
                parts = line.strip().split(':')
                if len(parts) == 2:
                    ele_id_str, nodes_str = parts
                    ele_id = int(ele_id_str.replace('Element', '').strip())
                    node_ids = [int(x) for x in nodes_str.replace('nodes', '').strip().split()]
                    if len(node_ids) == 4:
                        elements.append((ele_id, node_ids))
    except Exception as e:
        print(f"Error reading elements from {filename}: {e}")
    return elements

def parse_displacement_file(filename):
    """Parse the eq_displacements.out file to extract displacement data"""
    if not os.path.exists(filename):
        print(f"ERROR: {filename} not found.")
        return None, None
    try:
        data = np.loadtxt(filename)
        if data.size == 0:
            print("ERROR: Displacement file is empty.")
            return None, None
        time = data[:, 0]
        disp_data = data[:, 1:]
        return time, disp_data
    except Exception as e:
        print(f"Error reading displacement data from {filename}: {e}")
        return None, None

def parse_element_stress(stress_filename, elements):
    """Parse the eq_element_stress.out file to extract stress data for Q4 elements"""
    if not os.path.exists(stress_filename):
        print(f"ERROR: {stress_filename} not found.")
        return None
    try:
        data = np.loadtxt(stress_filename)
        if data.size == 0:
            print("ERROR: Stress file is empty.")
            return None
        
        # For Q4 elements, stress data format varies
        # Let's determine the actual format from the data
        if data.ndim == 1:
            # Single time step
            stress_data = data[1:].reshape(1, -1)  # Skip time column
        else:
            # Multiple time steps
            stress_data = data[:, 1:]  # Skip time column
        
        num_elements = len(elements)
        num_time_steps = stress_data.shape[0]
        stress_values_per_element = stress_data.shape[1] // num_elements
        
        print(f"Debug: {num_elements} elements, {num_time_steps} time steps, {stress_values_per_element} stress values per element")
        
        # Reshape to (time_steps, elements, stress_components)
        stress_reshaped = stress_data.reshape(num_time_steps, num_elements, stress_values_per_element)
        
        return stress_reshaped
    except Exception as e:
        print(f"Error reading stress data from {stress_filename}: {e}")
        return None

def plot_displacement_history():
    """Plot displacement history for the top node"""
    disp_file = 'RESULTS/eq_displacements.out'
    if not os.path.exists(disp_file):
        print("ERROR: Displacement file not found.")
        return
    time, disp_data = parse_displacement_file(disp_file)
    if time is None or disp_data is None:
        print("Skipping displacement plot due to missing or empty data.")
        return
    # Plot displacement for the top node (assuming it's the last node)
    top_node_idx = disp_data.shape[1] // 2 - 1  # Index for x-displacement of top node
    if top_node_idx >= 0 and 2*top_node_idx + 1 < disp_data.shape[1]:
        x_disp = disp_data[:, 2*top_node_idx]
        y_disp = disp_data[:, 2*top_node_idx + 1]
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 1, 1)
        plt.plot(time, x_disp, 'b-', linewidth=2)
        plt.title('Top Node Displacement History')
        plt.ylabel('X-Displacement (m)')
        plt.grid(True, alpha=0.3)
        plt.subplot(2, 1, 2)
        plt.plot(time, y_disp, 'r-', linewidth=2)
        plt.xlabel('Time (s)')
        plt.ylabel('Y-Displacement (m)')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('RESULTS/eq_displacement_history.png', dpi=300, bbox_inches='tight')
        plt.show()

def plot_mesh_and_deformed():
    """Plot original and deformed mesh side by side"""
    nodes_file = 'RESULTS/eq_nodes.txt'
    elements_file = 'RESULTS/eq_elements.txt'
    disp_file = 'RESULTS/eq_displacements.out'
    if not (os.path.exists(nodes_file) and os.path.exists(elements_file) and os.path.exists(disp_file)):
        print("ERROR: Required files for mesh plot not found.")
        return
    nodes = parse_nodes_file(nodes_file)
    elements = parse_elements_file(elements_file)
    time, disp_data = parse_displacement_file(disp_file)
    if time is None or disp_data is None:
        print("Skipping mesh plot due to missing or empty displacement data.")
        return
    # Get last time step displacements
    disp_last = disp_data[-1, :] if disp_data.shape[0] > 0 else np.zeros(disp_data.shape[1])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    # Original mesh
    node_dict = {node[0]: (node[1], node[2]) for node in nodes}
    x_coords = [node[1] for node in nodes]
    y_coords = [node[2] for node in nodes]
    ax1.scatter(x_coords, y_coords, c='black', s=50, zorder=5)
    for ele_id, node_ids in elements:
        if len(node_ids) == 4:
            coords = [node_dict[nid] for nid in node_ids if nid in node_dict]
            if len(coords) == 4:
                poly_x = [c[0] for c in coords]
                poly_y = [c[1] for c in coords]
                ax1.fill(poly_x, poly_y, alpha=0.3, color='lightblue', edgecolor='blue', linewidth=1)
                for i in range(4):
                    x1, y1 = coords[i]
                    x2, y2 = coords[(i + 1) % 4]
                    ax1.plot([x1, x2], [y1, y2], 'b-', linewidth=1, zorder=3)
    ax1.set_title('Original Mesh')
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')
    # Deformed mesh
    deformed_coords = {}
    for i, node in enumerate(nodes):
        if 2*i < len(disp_last):
            def_x = node[1] + disp_last[2*i]
            def_y = node[2] + disp_last[2*i + 1]
        else:
            def_x = node[1]
            def_y = node[2]
        deformed_coords[node[0]] = (def_x, def_y)
    def_x_coords = [deformed_coords[node[0]][0] for node in nodes]
    def_y_coords = [deformed_coords[node[0]][1] for node in nodes]
    ax2.scatter(def_x_coords, def_y_coords, c='black', s=50, zorder=5)
    for ele_id, node_ids in elements:
        if len(node_ids) == 4:
            coords = [deformed_coords[nid] for nid in node_ids if nid in deformed_coords]
            if len(coords) == 4:
                poly_x = [c[0] for c in coords]
                poly_y = [c[1] for c in coords]
                ax2.fill(poly_x, poly_y, alpha=0.3, color='lightcoral', edgecolor='red', linewidth=1)
                for i in range(4):
                    x1, y1 = coords[i]
                    x2, y2 = coords[(i + 1) % 4]
                    ax2.plot([x1, x2], [y1, y2], 'r-', linewidth=1, zorder=3)
    ax2.set_title('Deformed Mesh (last step)')
    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Y (m)')
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal')
    plt.tight_layout()
    plt.savefig('RESULTS/eq_mesh_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_stress_distribution():
    """Plot sigma_xx, sigma_yy, and tau_xy stress distributions on undeformed mesh as subplots"""
    elements = parse_elements_file('RESULTS/eq_elements.txt')
    if not elements:
        print("ERROR: No elements found in eq_elements.txt. Skipping stress plot.")
        return
    stress = parse_element_stress('RESULTS/eq_element_stress.out', elements)
    if stress is None:
        print("ERROR: No stress data found in eq_element_stress.out. The analysis may have failed.")
        return
    
    # Get the last time step stress data
    stress_last = stress[-1, :, :] if stress.shape[0] > 0 else stress[0, :, :]
    
    # Determine how many stress components we have
    num_stress_components = stress_last.shape[1]
    print(f"Debug: Found {num_stress_components} stress components per element")
    
    # Use only the first 3 components (sigma_xx, sigma_yy, tau_xy) if we have more
    if num_stress_components >= 3:
        stress_components = 3
        stress_names = [r'$\sigma_{xx}$', r'$\sigma_{yy}$', r'$\tau_{xy}$']
    else:
        stress_components = num_stress_components
        stress_names = [f'Stress_{i+1}' for i in range(stress_components)]
    
    nodes = parse_nodes_file('RESULTS/eq_nodes.txt')
    node_dict = {node[0]: (node[1], node[2]) for node in nodes}
    all_x = [node[1] for node in nodes]
    all_y = [node[2] for node in nodes]
    xlim = (min(all_x), max(all_x))
    ylim = (min(all_y), max(all_y))
    
    fig, axes = plt.subplots(stress_components, 1, figsize=(10, 8), sharex=True, sharey=True)
    if stress_components == 1:
        axes = [axes]
    
    for k, ax in enumerate(axes):
        cmap = plt.get_cmap('jet')
        vals = stress_last[:, k]
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
                    centroid = (np.mean(poly_x), np.mean(poly_y))
                    centroids.append(centroid)
                else:
                    centroids.append((None, None))
            else:
                centroids.append((None, None))
        ax.set_title(f'{stress_names[k]} Stress')
        ax.set_xlabel('X (m)')
        if k == 0:
            ax.set_ylabel('Y (m)')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.1)
        cb = fig.colorbar(sm, cax=cax, label=stress_names[k] + ' (Pa)')
        if centroids[max_idx][0] is not None:
            ax.annotate(f"max: {vmax:.2e}", xy=centroids[max_idx], xytext=(centroids[max_idx][0], centroids[max_idx][1]),
                        color='red', fontsize=10, fontweight='bold', ha='center', va='center',
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="red", lw=1, alpha=0.7))
        if centroids[min_idx][0] is not None:
            ax.annotate(f"min: {vmin:.2e}", xy=centroids[min_idx], xytext=(centroids[min_idx][0], centroids[min_idx][1]),
                        color='blue', fontsize=10, fontweight='bold', ha='center', va='center',
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="blue", lw=1, alpha=0.7))
    plt.tight_layout()
    plt.savefig('RESULTS/eq_stress_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()

def animate_deformation():
    """Generate and save an animation of mesh deformation over time."""
    disp_file = 'RESULTS/eq_displacements.out'
    nodes_file = 'RESULTS/eq_nodes.txt'
    elements_file = 'RESULTS/eq_elements.txt'
    if not (os.path.exists(disp_file) and os.path.exists(nodes_file) and os.path.exists(elements_file)):
        print("ERROR: Required files for animation not found.")
        return
    # Parse nodes and elements
    nodes = parse_nodes_file(nodes_file)
    elements = parse_elements_file(elements_file)
    time, disp_data = parse_displacement_file(disp_file)
    if time is None or disp_data is None:
        print("Skipping animation due to missing or empty displacement data.")
        return
    # Prepare node order for fast lookup
    node_id_to_idx = {node[0]: idx for idx, node in enumerate(nodes)}
    node_coords = np.array([[node[1], node[2]] for node in nodes])
    n_frames = disp_data.shape[0]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_aspect('equal')
    # Get mesh limits
    x_coords = node_coords[:, 0]
    y_coords = node_coords[:, 1]
    ax.set_xlim(np.min(x_coords), np.max(x_coords))
    ax.set_ylim(np.min(y_coords), np.max(y_coords))
    # Plot undeformed mesh (light gray)
    for ele_id, node_ids in elements:
        if len(node_ids) == 4:
            coords = [node_coords[node_id_to_idx[nid]] for nid in node_ids]
            poly_x = [c[0] for c in coords]
            poly_y = [c[1] for c in coords]
            ax.fill(poly_x, poly_y, alpha=0.1, color='gray', edgecolor='gray', linewidth=0.5)
    # Deformed mesh (animated)
    deformed_patches = []
    for ele_id, node_ids in elements:
        if len(node_ids) == 4:
            coords = [node_coords[node_id_to_idx[nid]] for nid in node_ids]
            patch = plt.Polygon(coords, closed=True, edgecolor='blue', facecolor='lightblue', alpha=0.7, lw=1)
            ax.add_patch(patch)
            deformed_patches.append((patch, node_ids))
    time_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, va='top', ha='left', fontsize=12, color='black')
    def update(frame):
        disp = disp_data[frame, :]
        for patch, node_ids in deformed_patches:
            coords = []
            for nid in node_ids:
                idx = node_id_to_idx[nid]
                x = node_coords[idx, 0] + disp[2*idx]
                y = node_coords[idx, 1] + disp[2*idx+1]
                coords.append([x, y])
            patch.set_xy(coords)
        time_text.set_text(f"Time: {time[frame]:.3f} s")
        return [p[0] for p in deformed_patches] + [time_text]
    ani = animation.FuncAnimation(fig, update, frames=n_frames, interval=50, blit=True)
    output_path = 'RESULTS/eq_deformation_animation.mp4'
    ani.save(output_path, writer='ffmpeg', fps=20)
    print(f"Animation saved as {output_path}")
    plt.close(fig)

def main():
    print("=== OpenSees Earthquake Model Results Visualization ===")
    if not os.path.exists('RESULTS'):
        print("ERROR: RESULTS directory not found. Run the OpenSees analysis first.")
        return
    print("Creating top node displacement history plot...")
    plot_displacement_history()
    print("Creating mesh comparison plot...")
    plot_mesh_and_deformed()
    print("Creating stress distribution plot...")
    plot_stress_distribution()
    print("Creating deformation animation...")
    animate_deformation()
    print("\nVisualization complete!")
    print("Check the RESULTS directory for output images and animation:")
    print("  - eq_displacement_history.png")
    print("  - eq_mesh_comparison.png")
    print("  - eq_stress_distribution.png")
    print("  - eq_deformation_animation.mp4")

if __name__ == "__main__":
    main() 