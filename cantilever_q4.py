#!/usr/bin/env python3
"""
2D Q4 Element Cantilever Beam Analysis
Created for simplified OpenSees (Linux serial build)
"""

import os
import sys
import numpy as np

# Import OpenSees
import openseespy.opensees as ops


def main():
    """Main function to perform cantilever beam analysis"""
    
    # Clear any previous model
    ops.wipe()
    
    # Create the model
    ops.model('BasicBuilder', '-ndm', 2, '-ndf', 2)
    
    # Material properties
    E = 30000.0    # Young's modulus (ksi)
    nu = 0.3       # Poisson's ratio
    rho = 0.0      # Density (negligible for static analysis)
    
    # Create elastic material
    ops.nDMaterial('ElasticIsotropic', 1, E, nu, rho)
    
    # Beam geometry
    L = 100.0      # Length (in)
    H = 10.0       # Height (in)
    t = 1.0        # Thickness (in)
    
    # Mesh parameters
    nx = 20        # Number of elements in x-direction
    ny = 4         # Number of elements in y-direction
    dx = L / nx
    dy = H / ny
    
    print("=== Cantilever Beam Model Information ===")
    print(f"Beam dimensions: {L} x {H} x {t} inches")
    print(f"Mesh: {nx} x {ny} elements ({dx:.2f} x {dy:.2f} inch elements)")
    print(f"Material: E = {E} ksi, nu = {nu}, rho = {rho}")
    print(f"Total nodes: {(nx+1) * (ny+1)}")
    print(f"Total elements: {nx * ny}")
    print("=" * 40)
    
    # Create nodes
    nodeID = 1
    nodes = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x = i * dx
            y = j * dy
            ops.node(nodeID, x, y)
            nodes.append((nodeID, x, y))
            nodeID += 1
    
    # Create Q4 elements
    eleID = 1
    elements = []
    for j in range(ny):
        for i in range(nx):
            n1 = i + j * (nx + 1) + 1
            n2 = i + 1 + j * (nx + 1) + 1
            n3 = i + 1 + (j + 1) * (nx + 1) + 1
            n4 = i + (j + 1) * (nx + 1) + 1
            
            ops.element('quad', eleID, n1, n2, n3, n4, t, 'PlaneStress', 1)
            elements.append((eleID, n1, n2, n3, n4))
            eleID += 1
    
    # Apply boundary conditions (fixed at left end)
    fixed_nodes = []
    for j in range(ny + 1):
        nodeID = j * (nx + 1) + 1
        ops.fix(nodeID, 1, 1)
        fixed_nodes.append(nodeID)
    
    print(f"Fixed nodes (left end): {fixed_nodes}")
    
    # Create load pattern
    ops.pattern('Plain', 1, 'Linear')
    
    # Apply tip load at right end (distributed among top nodes)
    loadPerNode = -100.0 / (ny + 1)
    loaded_nodes = []
    for j in range(ny + 1):
        nodeID = j * (nx + 1) + nx + 1
        ops.load(nodeID, 0.0, loadPerNode)
        loaded_nodes.append(nodeID)
    
    print(f"Loaded nodes (right end): {loaded_nodes}")
    print(f"Load per node: {loadPerNode:.2f} kips")
    
    # Create analysis
    ops.system('UmfPack')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.numberer('RCM')
    ops.analysis('Static')
    
    # Create RESULTS directory
    results_dir = 'RESULTS'
    os.makedirs(results_dir, exist_ok=True)
    
    # Create recorders
    ops.recorder('Node', f'{results_dir}/displacements.out', 'disp', '-time', '-nodes', 'all')
    ops.recorder('Element', f'{results_dir}/elementForces.out', 'force', '-time', '-ele', 'all')
    ops.recorder('Element', f'{results_dir}/elementStrains.out', 'material', '-time', '-ele', 'all')
    ops.recorder('Node', f'{results_dir}/reactions.out', 'reaction', '-time', '-nodes', 'all')

    # Save model information similar to TCL print commands
    ops.printModel('-file', f'{results_dir}/nodes.txt', '-node')
    ops.printModel('-file', f'{results_dir}/elements.txt', '-ele')
    ops.printModel('-file', f'{results_dir}/materials.txt', '-material')
    ops.printModel('-file', f'{results_dir}/patterns.txt', '-pattern')
    ops.printModel('-file', f'{results_dir}/constraints.txt', '-constraint')
    
    # Perform analysis
    print("\n=== Starting Analysis ===")
    result = ops.analyze(1)
    
    if result == 0:
        print("Analysis completed successfully!")
    else:
        print(f"Analysis failed with error code: {result}")
        return
    
    # Print results summary
    print("\n=== Results Summary ===")
    print(f"Results saved in {results_dir} directory:")
    print("  - displacements.out: Node displacements")
    print("  - elementForces.out: Element forces")
    print("  - elementStrains.out: Element strains")
    print("  - reactions.out: Support reactions")
    
    # Print tip displacement
    tipNode = ny * (nx + 1) + nx + 1
    tipDisp = ops.nodeDisp(tipNode, 2)
    print(f"Tip displacement: {tipDisp:.6f} inches")
    
    # Print model statistics
    print("\n=== Model Statistics ===")
    print(f"Total nodes created: {len(nodes)}")
    print(f"Total elements created: {len(elements)}")
    print(f"Fixed nodes: {len(fixed_nodes)}")
    print(f"Loaded nodes: {len(loaded_nodes)}")
    
    # Save model information to file
    model_info_file = f'{results_dir}/model_info.txt'
    with open(model_info_file, 'w') as f:
        f.write("=== Cantilever Beam Model Information ===\n")
        f.write(f"Beam dimensions: {L} x {H} x {t} inches\n")
        f.write(f"Mesh: {nx} x {ny} elements ({dx:.2f} x {dy:.2f} inch elements)\n")
        f.write(f"Material: E = {E} ksi, nu = {nu}, rho = {rho}\n")
        f.write(f"Total nodes: {len(nodes)}\n")
        f.write(f"Total elements: {len(elements)}\n")
        f.write(f"Fixed nodes: {fixed_nodes}\n")
        f.write(f"Loaded nodes: {loaded_nodes}\n")
        f.write(f"Load per node: {loadPerNode:.2f} kips\n")
        f.write(f"Tip displacement: {tipDisp:.6f} inches\n")
    
    print(f"Model information saved to: {model_info_file}")
    
    # Clean up
    ops.wipe()
    print("\nModel cleared from memory.")

if __name__ == "__main__":
    main() 