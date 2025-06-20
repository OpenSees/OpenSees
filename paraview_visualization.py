#!/usr/bin/env python3
"""
ParaView Visualization Script for Cantilever Beam Analysis
This script creates a ParaView visualization of the OpenSees results
"""

import numpy as np
import matplotlib.pyplot as plt
import os

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

def create_paraview_visualization():
    """Create ParaView visualization of cantilever beam results"""
    
    # Check if RESULTS directory exists
    if not os.path.exists('RESULTS'):
        print("ERROR: RESULTS directory not found. Run the OpenSees analysis first.")
        return
    
    try:
        from paraview import simple
    except ImportError:
        print("ParaView not available, using matplotlib instead...")
        create_matplotlib_visualization()
        return
    
    # Create ParaView pipeline
    print("Creating ParaView visualization...")
    
    # Create a new render view
    renderView = simple.GetActiveViewOrCreate('RenderView')
    renderView.ViewSize = [1200, 800]
    renderView.CenterOfRotation = [50.0, 5.0, 0.0]
    renderView.CameraPosition = [50.0, 5.0, 100.0]
    renderView.CameraFocalPoint = [50.0, 5.0, 0.0]
    renderView.CameraViewUp = [0.0, 1.0, 0.0]
    
    # Parse displacement data
    time_step, disp_data = parse_displacement_file('RESULTS/displacements.out')
    nodes = parse_nodes_file('RESULTS/nodes.txt')
    
    # Create a CSV file for ParaView
    csv_filename = 'RESULTS/paraview_data.csv'
    with open(csv_filename, 'w') as f:
        f.write("X,Y,Z,Displacement_X,Displacement_Y,Displacement_Magnitude\n")
        for i, node in enumerate(nodes):
            if 2*i < len(disp_data):
                x = node[1]
                y = node[2]
                z = 0.0
                disp_x = disp_data[2*i]
                disp_y = disp_data[2*i + 1]
                disp_mag = np.sqrt(disp_x**2 + disp_y**2)
                f.write(f"{x},{y},{z},{disp_x},{disp_y},{disp_mag}\n")
    
    # Load the CSV data
    csvReader = simple.CSVReader(FileName=[csv_filename])
    csvReader.HaveHeaders = True
    
    # Create a table to points filter
    tableToPoints = simple.TableToPoints(Input=csvReader)
    tableToPoints.XColumn = 'X'
    tableToPoints.YColumn = 'Y'
    tableToPoints.ZColumn = 'Z'
    
    # Create a glyph3D to show displacement vectors
    glyph = simple.Glyph3D(Input=tableToPoints)
    glyph.GlyphType = 'Arrow'
    glyph.ScaleFactor = 5.0  # Scale factor for displacement visualization
    glyph.GlyphMode = 'All Points'
    glyph.Vectors = ['POINTS', 'Displacement_X', 'Displacement_Y', 'Displacement_Magnitude']
    
    # Display the glyph
    glyphDisplay = simple.Show(glyph, renderView)
    glyphDisplay.Representation = 'Surface'
    glyphDisplay.ColorArrayName = 'Displacement_Magnitude'
    glyphDisplay.LookupTable = simple.GetColorTransferFunction('Displacement_Magnitude')
    glyphDisplay.SetScalarBarVisibility(renderView, True)
    
    # Set up color map
    colorMap = simple.GetColorTransferFunction('Displacement_Magnitude')
    colorMap.RGBPoints = [0.0, 0.0, 0.0, 1.0, 0.5, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0]
    colorMap.ColorSpace = 'RGB'
    
    # Add text annotation for title
    text = simple.Text()
    text.Text = 'Cantilever Beam Analysis Results'
    textDisplay = simple.Show(text, renderView)
    textDisplay.Position = [0.1, 0.9]
    textDisplay.FontSize = 16
    
    # Reset camera to fit all data
    simple.ResetCamera()
    
    # Save screenshot
    simple.SaveScreenshot('RESULTS/beam_visualization.png', renderView)
    
    print("Visualization created successfully!")
    print("Screenshot saved as: RESULTS/beam_visualization.png")
    print("Use matplotlib visualization for detailed analysis of results")

def create_matplotlib_visualization():
    """Create matplotlib visualization as an alternative to ParaView"""
    
    print("Creating matplotlib visualization...")
    
    # Parse displacement data
    time_step, disp_data = parse_displacement_file('RESULTS/displacements.out')
    nodes = parse_nodes_file('RESULTS/nodes.txt')
    
    if len(nodes) == 0:
        print("ERROR: No nodes found in nodes.txt")
        return
    
    # Parse elements if available
    elements = []
    if os.path.exists('RESULTS/elements.txt'):
        elements = parse_elements_file('RESULTS/elements.txt')
    
    # Extract coordinates and displacements
    x_coords = np.array([node[1] for node in nodes])
    y_coords = np.array([node[2] for node in nodes])
    
    # Extract displacements
    disp_x = []
    disp_y = []
    for i in range(0, len(disp_data), 2):
        if i + 1 < len(disp_data):
            disp_x.append(disp_data[i])
            disp_y.append(disp_data[i + 1])
    
    # Ensure we have the same number of displacement values as nodes
    if len(disp_x) != len(nodes):
        print(f"Warning: Number of displacement values ({len(disp_x)}) doesn't match number of nodes ({len(nodes)})")
        # Pad or truncate to match
        if len(disp_x) < len(nodes):
            disp_x.extend([0.0] * (len(nodes) - len(disp_x)))
            disp_y.extend([0.0] * (len(nodes) - len(disp_y)))
        else:
            disp_x = disp_x[:len(nodes)]
            disp_y = disp_y[:len(nodes)]
    
    disp_x = np.array(disp_x)
    disp_y = np.array(disp_y)
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Create node dictionary for easy lookup
    node_dict = {node[0]: (node[1], node[2]) for node in nodes}
    
    # Plot 1: Original mesh
    # Plot nodes
    ax1.scatter(x_coords, y_coords, c='black', s=50, zorder=5)
    
    # Plot elements as filled polygons if available
    if elements:
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
    
    # Plot 2: Deformed mesh
    # Calculate deformed coordinates
    deformed_coords = {}
    for i, node in enumerate(nodes):
        if i < len(disp_x):
            def_x = node[1] + disp_x[i]
            def_y = node[2] + disp_y[i]
        else:
            def_x = node[1]
            def_y = node[2]
        deformed_coords[node[0]] = (def_x, def_y)
    
    # Plot deformed nodes
    def_x_coords = [deformed_coords[node[0]][0] for node in nodes]
    def_y_coords = [deformed_coords[node[0]][1] for node in nodes]
    ax2.scatter(def_x_coords, def_y_coords, c='black', s=50, zorder=5)
    
    # Plot deformed elements if available
    if elements:
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
    plt.savefig('RESULTS/beam_matplotlib.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Matplotlib visualization saved as: RESULTS/beam_matplotlib.png")

def main():
    """Main function"""
    print("=== ParaView Visualization for Cantilever Beam ===")
    
    # Try ParaView visualization first
    try:
        create_paraview_visualization()
    except ImportError:
        print("ParaView not available, using matplotlib instead...")
        create_matplotlib_visualization()
    except Exception as e:
        print(f"ParaView visualization failed: {e}")
        print("Falling back to matplotlib visualization...")
        create_matplotlib_visualization()
    
    print("\nVisualization complete!")
    print("Check the RESULTS directory for output files.")

if __name__ == "__main__":
    main() 