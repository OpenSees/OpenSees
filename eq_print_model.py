import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def draw_dashpot(ax, p1, p2, size=0.1, piston_width=0.05):
    """
    Draws a dashpot symbol between two points p1 and p2 on the given axis.

    Args:
        ax: The matplotlib axis to draw on.
        p1 (tuple): The (x, y) coordinates of the starting point.
        p2 (tuple): The (x, y) coordinates of the ending point.
        size (float): The length of the dashpot symbol's main body.
        piston_width (float): The width of the piston head.
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    
    # Vector and direction
    vec = p2 - p1
    dist = np.linalg.norm(vec)
    if dist == 0: return
    direction = vec / dist
    perp_dir = np.array([-direction[1], direction[0]]) # Perpendicular vector

    # Dashpot body (cylinder)
    body_start = p1 + direction * (dist * 0.5 - size)
    body_end = p1 + direction * (dist * 0.5 + size)
    
    # Cylinder lines
    p_bl = body_start - perp_dir * piston_width
    p_br = body_start + perp_dir * piston_width
    p_tl = body_end - perp_dir * piston_width
    p_tr = body_end + perp_dir * piston_width

    ax.plot([p_bl[0], p_tl[0]], [p_bl[1], p_tl[1]], color='gray', linewidth=1)
    ax.plot([p_br[0], p_tr[0]], [p_br[1], p_tr[1]], color='gray', linewidth=1)
    ax.plot([p_tl[0], p_tr[0]], [p_tl[1], p_tr[1]], color='gray', linewidth=1)

    # Piston line and head
    piston_start = p1 + direction * (dist * 0.1)
    piston_end = p1 + direction * (dist * 0.5)
    ax.plot([piston_start[0], piston_end[0]], [piston_start[1], piston_end[1]], color='gray', linewidth=1)
    
    head_left = piston_end - perp_dir * (piston_width * 1.5)
    head_right = piston_end + perp_dir * (piston_width * 1.5)
    ax.plot([head_left[0], head_right[0]], [head_left[1], head_right[1]], color='gray', linewidth=1.5)

    # Connection lines
    ax.plot([p1[0], piston_start[0]], [p1[1], piston_start[1]], color='gray', linewidth=1, linestyle='--')
    ax.plot([body_start[0], p2[0]], [body_start[1], p2[1]], color='gray', linewidth=1, linestyle='--')

def plot_model_from_parameters():
    # Parameters (should match earthquake_model.tcl)
    nx = 10
    ny = 5
    width = 10.0
    height = 5.0
    thickness = 1.0
    dx = width / nx
    dy = height / ny
    # Generate nodes
    nodes = {}
    node_id = 1
    for j in range(ny+1):
        for i in range(nx+1):
            x = i * dx
            y = j * dy
            nodes[node_id] = (x, y)
            node_id += 1
    # Generate quad elements
    quad_elements = []
    for j in range(ny):
        for i in range(nx):
            n1 = i + j*(nx+1) + 1
            n2 = i+1 + j*(nx+1) + 1
            n3 = i+1 + (j+1)*(nx+1) + 1
            n4 = i + (j+1)*(nx+1) + 1
            quad_elements.append([n1, n2, n3, n4])
    # Far-field node generation (for dashpots)
    far_field_nodes = {}
    dashpot_pairs = []
    ff_id = (nx+1)*(ny+1) + 1
    # LEFT and RIGHT boundaries
    for j in range(ny+1):
        # Left
        nodeL = j*(nx+1) + 1
        far_field_nodes[ff_id] = nodes[nodeL]
        dashpot_pairs.append((nodeL, ff_id, 'x'))
        dashpot_pairs.append((nodeL, ff_id, 'y'))
        ff_id += 1
        # Right
        nodeR = j*(nx+1) + nx + 1
        far_field_nodes[ff_id] = nodes[nodeR]
        dashpot_pairs.append((nodeR, ff_id, 'x'))
        dashpot_pairs.append((nodeR, ff_id, 'y'))
        ff_id += 1
    # BOTTOM boundary (excluding corners)
    for i in range(1, nx):
        nodeB = i + 1
        far_field_nodes[ff_id] = nodes[nodeB]
        dashpot_pairs.append((nodeB, ff_id, 'x'))
        dashpot_pairs.append((nodeB, ff_id, 'y'))
        ff_id += 1
    # Plotting
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_aspect('equal')
    # Plot quad elements
    for ele_nodes in quad_elements:
        poly_coords = [nodes[nid] for nid in ele_nodes]
        polygon = patches.Polygon(poly_coords, closed=True, edgecolor='black', facecolor='skyblue', alpha=0.7, lw=1)
        ax.add_patch(polygon)
    # Plot dashpot elements
    for n1, n2, direction in dashpot_pairs:
        p1 = nodes[n1]
        p2 = far_field_nodes[n2]
        draw_dashpot(ax, p1, p2, size=min(dx, dy)*0.4, piston_width=0.1)
    # Plot mesh nodes
    x_coords, y_coords = zip(*nodes.values())
    ax.plot(x_coords, y_coords, 'ko', markersize=3, label='Mesh Nodes')
    # Plot far-field nodes
    if far_field_nodes:
        ff_x, ff_y = zip(*far_field_nodes.values())
        ax.plot(ff_x, ff_y, 'rs', markersize=5, label='Fixed Far-Field Nodes')
    # Plot earthquake load arrow at a single top surface node (center of top boundary)
    center_i = nx // 2
    top_node = ny * (nx + 1) + center_i + 1
    x, y = nodes[top_node]
    arrow_length = height * 0.12
    ax.arrow(x, y, arrow_length, 0, head_width=0.12*dy, head_length=0.18*arrow_length, fc='red', ec='red', linewidth=2)
    # Add a single legend entry for earthquake load
    ax.plot([], [], color='red', linewidth=2, label='Earthquake Load (Shear, single top node)')
    ax.set_xlabel("X-coordinate (m)")
    ax.set_ylabel("Y-coordinate (m)")
    ax.set_title("OpenSees Model: Mesh and Half-Space Boundaries (from parameters)")
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend()
    plt.tight_layout()
    plt.savefig('RESULTS/eq_model.png', dpi=300, bbox_inches='tight')
    print("Model plot saved as RESULTS/eq_model.png")
    # plt.show()  # Optionally, comment out to avoid showing interactively

if __name__ == "__main__":
    plot_model_from_parameters()
