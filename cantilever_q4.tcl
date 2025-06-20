# 2D Q4 Element Cantilever Beam Analysis
# Created for simplified OpenSees (Linux serial build)

# Clear any previous model
wipe

# Create the model
model BasicBuilder -ndm 2 -ndf 2

# Material properties
set E 30000.0    ;# Young's modulus (ksi)
set nu 0.3       ;# Poisson's ratio
set rho 0.0      ;# Density (negligible for static analysis)

# Create elastic material
nDMaterial ElasticIsotropic 1 $E $nu $rho

# Beam geometry
set L 100.0      ;# Length (in)
set H 10.0       ;# Height (in)
set t 1.0        ;# Thickness (in)

# Mesh parameters
set nx 5        ;# Number of elements in x-direction
set ny 2         ;# Number of elements in y-direction
set dx [expr $L/$nx]
set dy [expr $H/$ny]

# Create nodes
set nodeID 1
for {set j 0} {$j <= $ny} {incr j} {
    for {set i 0} {$i <= $nx} {incr i} {
        set x [expr $i * $dx]
        set y [expr $j * $dy]
        node $nodeID $x $y
        incr nodeID
    }
}

# Create Q4 elements
set eleID 1
for {set j 0} {$j < $ny} {incr j} {
    for {set i 0} {$i < $nx} {incr i} {
        set n1 [expr $i + $j*($nx+1) + 1]
        set n2 [expr $i + 1 + $j*($nx+1) + 1]
        set n3 [expr $i + 1 + ($j+1)*($nx+1) + 1]
        set n4 [expr $i + ($j+1)*($nx+1) + 1]
        
        element quad $eleID $n1 $n2 $n3 $n4 $t "PlaneStress" 1
        incr eleID
    }
}

# Apply boundary conditions (fixed at left end)
for {set j 0} {$j <= $ny} {incr j} {
    set nodeID [expr $j*($nx+1) + 1]
    fix $nodeID 1 1
}

# Create load pattern
pattern Plain 1 "Linear" {
    # Apply tip load at right end (distributed among top nodes)
    set loadPerNode [expr -100.0 / ($ny + 1)]
    for {set j 0} {$j <= $ny} {incr j} {
        set nodeID [expr $j*($nx+1) + $nx + 1]
        load $nodeID 0.0 $loadPerNode
    }
}

# Create analysis
system UmfPack
constraints Plain
integrator LoadControl 1.0
algorithm Linear
numberer RCM
analysis Static

# Create RESULTS directory
file mkdir RESULTS

# Calculate total nodes and elements
set totalNodes [expr ($nx+1)*($ny+1)]
set totalElements [expr $nx * $ny]

# Print model information to files (simplified)
print -file RESULTS/model_info.txt
print -file RESULTS/nodes.txt -node
print -file RESULTS/elements.txt -ele

# Create recorders (simplified to match working version)
recorder Node -file RESULTS/displacements.out -time -nodeRange 1 $totalNodes -dof 1 2 disp

# Perform analysis
analyze 1

# Print results
puts "Analysis completed successfully!"
puts "Results saved in RESULTS directory:"
puts "  - model_info.txt: Complete model information"
puts "  - nodes.txt: Node coordinates and DOFs"
puts "  - elements.txt: Element connectivity and properties"
puts "  - displacements.out: Node displacements"

# Print tip displacement
set tipNode [expr $ny*($nx+1) + $nx + 1]
set tipDisp [nodeDisp $tipNode 2]
puts "Tip displacement: $tipDisp inches"

# Clean up
wipe 