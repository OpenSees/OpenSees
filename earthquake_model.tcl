# earthquake_model.tcl
# 2D Elastic Soil Model with Elastic Half-Space Boundary Conditions
# Compatible with OpenSees. Outputs results in RESULTS/ with 'eq_' prefix.
#
# Revision Notes:
# - Model updated to use regular Q4 elements (simplified from SSPquad)
# - Added viscous dashpots (absorbing boundaries) on the left, right, and
#   bottom of the mesh to simulate an infinite domain and prevent wave reflections.
# - Loading changed from UniformExcitation to a compliant base input, where
#   the earthquake is applied as a shear stress time history at the base nodes.

# ----------------------
# Model Initialization
# ----------------------
wipe
model BasicBuilder -ndm 2 -ndf 2

# ----------------------
# User Parameters
# ----------------------
# Soil and mesh parameters
set E_soil 3.0e7         ;# Young's Modulus (Pa)
set nu_soil 0.3          ;# Poisson's Ratio
set rho_soil 2000.0      ;# Density (kg/m^3)
set thickness 1.0        ;# Out-of-plane thickness (m)
set width 10.0           ;# Model width (m)
set height 5.0           ;# Model height (m)
set nx 10                ;# Number of elements in x-direction
set ny 5                 ;# Number of elements in y-direction

# Derived mesh spacing
set dx [expr $width / $nx]
set dy [expr $height / $ny]

# ----------------------
# Material Definition
# ----------------------
set soilMatTag 1
nDMaterial ElasticIsotropic $soilMatTag $E_soil $nu_soil $rho_soil

# ----------------------
# Node and Element Generation
# ----------------------
# Node Generation
set nodeID 1
for {set j 0} {$j <= $ny} {incr j} {
    for {set i 0} {$i <= $nx} {incr i} {
        node $nodeID [expr $i * $dx] [expr $j * $dy]
        incr nodeID
    }
}

# Element Generation (using regular Q4 elements)
set meshElementList {}
set eleID 1
for {set j 0} {$j < $ny} {incr j} {
    for {set i 0} {$i < $nx} {incr i} {
        set n1 [expr $i + $j*($nx+1) + 1]
        set n2 [expr $i + 1 + $j*($nx+1) + 1]
        set n3 [expr $i + 1 + ($j+1)*($nx+1) + 1]
        set n4 [expr $i + ($j+1)*($nx+1) + 1]
        element quad $eleID $n1 $n2 $n3 $n4 $thickness PlaneStrain $soilMatTag
        lappend meshElementList $eleID
        incr eleID
    }
}

# ----------------------------------------------------
# Boundary Conditions - Elastic Half-Space
# ----------------------------------------------------
# This section implements viscous dashpots on the boundaries to absorb
# outgoing waves, simulating an infinite domain.

# Calculate wave velocities
set G_soil [expr {$E_soil / (2.0 * (1.0 + $nu_soil))}]
set Vs_soil [expr {sqrt($G_soil / $rho_soil)}]
set Vp_soil [expr {$Vs_soil * sqrt(2.0*(1.0-$nu_soil)/(1.0-2.0*$nu_soil))}]

# Start far-field node IDs after all mesh nodes
set farFieldNodeID [expr ($nx+1)*($ny+1) + 1]
set dashpotID 1000
set matID 5000

# LEFT and RIGHT boundaries
for {set j 0} {$j <= $ny} {incr j} {
    set tributaryArea [expr $dy * $thickness]
    if {$j == 0 || $j == $ny} { set tributaryArea [expr $tributaryArea / 2.0] }

    # Left boundary node
    set nodeL [expr $j*($nx+1) + 1]
    node $farFieldNodeID [nodeCoord $nodeL 1] [nodeCoord $nodeL 2]
    fix $farFieldNodeID 1 1
    # P-wave dashpot (normal, x-dir)
    set matP $matID
    uniaxialMaterial Viscous $matP [expr $rho_soil * $Vp_soil * $tributaryArea] 1.0
    element zeroLength $dashpotID $nodeL $farFieldNodeID -mat $matP -dir 1
    incr matID
    incr dashpotID
    # S-wave dashpot (tangential, y-dir)
    set matS $matID
    uniaxialMaterial Viscous $matS [expr $rho_soil * $Vs_soil * $tributaryArea] 1.0
    element zeroLength $dashpotID $nodeL $farFieldNodeID -mat $matS -dir 2
    incr matID
    incr dashpotID
    incr farFieldNodeID

    # Right boundary node
    set nodeR [expr $j*($nx+1) + $nx + 1]
    node $farFieldNodeID [nodeCoord $nodeR 1] [nodeCoord $nodeR 2]
    fix $farFieldNodeID 1 1
    set matP $matID
    uniaxialMaterial Viscous $matP [expr $rho_soil * $Vp_soil * $tributaryArea] 1.0
    element zeroLength $dashpotID $nodeR $farFieldNodeID -mat $matP -dir 1
    incr matID
    incr dashpotID
    set matS $matID
    uniaxialMaterial Viscous $matS [expr $rho_soil * $Vs_soil * $tributaryArea] 1.0
    element zeroLength $dashpotID $nodeR $farFieldNodeID -mat $matS -dir 2
    incr matID
    incr dashpotID
    incr farFieldNodeID
}

# BOTTOM boundary (excluding corners already done)
for {set i 1} {$i < $nx} {incr i} {
    set tributaryArea [expr $dx * $thickness]
    set nodeB [expr $i + 1]
    node $farFieldNodeID [nodeCoord $nodeB 1] [nodeCoord $nodeB 2]
    fix $farFieldNodeID 1 1
    # S-wave dashpot (x-dir)
    set matS $matID
    uniaxialMaterial Viscous $matS [expr $rho_soil * $Vs_soil * $tributaryArea] 1.0
    element zeroLength $dashpotID $nodeB $farFieldNodeID -mat $matS -dir 1
    incr matID
    incr dashpotID
    # P-wave dashpot (y-dir)
    set matP $matID
    uniaxialMaterial Viscous $matP [expr $rho_soil * $Vp_soil * $tributaryArea] 1.0
    element zeroLength $dashpotID $nodeB $farFieldNodeID -mat $matP -dir 2
    incr matID
    incr dashpotID
    incr farFieldNodeID
}

# ----------------------------------------------------
# Earthquake Loading - Applied as Shear Stress
# ----------------------------------------------------
# For a compliant base, the input motion is applied as a shear stress
# time history at the base nodes. This simulates waves entering the model from below.

set motionTag 1
set tStart 0.0
set tEnd 5.0
set period 0.5
set accel_factor 50.0 ;# This is the amplitude of the input ACCELERATION (increased from 1.0 to 50.0)

# The input acceleration is a sine wave: a(t) = accel_factor * sin(w*t)
# The corresponding velocity is: v(t) = -(accel_factor/w) * cos(w*t)
set omega [expr {2.0 * 3.141592653589793 / $period}]
set vel_factor [expr {$accel_factor / $omega}]

# Create a velocity time series (cosine wave)
# We use a negative factor because the integral of sin(wt) is -(1/w)cos(wt)
timeSeries Trig $motionTag $tStart $tEnd $period -factor [expr -$vel_factor] -shift [expr {3.141592653589793 / 2.0}]

# Define a Plain pattern with the velocity time series
pattern Plain 1 $motionTag {
    # Apply shear forces to the base nodes. The force is F(t) = 2*rho*Vs*v(t)*Area
    # The '2' is because this is an "outcrop" motion applied to a compliant base.
    # The load command takes a load factor, which is multiplied by the time series.
    # So the load factor is 2 * rho * Vs * Area.
    for {set i 0} {$i <= $nx} {incr i} {
        set baseNodeID [expr $i + 1]
        set tributaryArea [expr $dx * $thickness]
        if {$i == 0 || $i == $nx} { set tributaryArea [expr $tributaryArea / 2.0] }

        set shearForceFactor [expr {2.0 * $rho_soil * $Vs_soil * $tributaryArea}]
        
        # Load is applied in the horizontal (DOF 1) direction
        load $baseNodeID $shearForceFactor 0.0
    }
}

# ----------------------
# Analysis Setup
# ----------------------
constraints Plain
numberer RCM
system UmfPack
integrator Newmark 0.5 0.25
algorithm Linear

# Rayleigh damping for internal material damping
set zeta 0.02
set f1 [expr {$Vs_soil/(4.0*$height)}]
set w1 [expr {2.0*3.141592653589793*$f1}]
set alphaM [expr {2.0*$zeta*$w1}]
set betaK 0.0
rayleigh $alphaM $betaK 0.0 0.0

analysis Transient

# ----------------------
# Output Setup
# ----------------------
file mkdir RESULTS
set totalNodes [expr ($nx+1)*($ny+1)]
set totalElements [expr $nx * $ny]

# Record displacements and accelerations for all mesh nodes
recorder Node -file RESULTS/eq_displacements.out -time -nodeRange 1 $totalNodes -dof 1 2 disp
recorder Node -file RESULTS/eq_accelerations.out -time -nodeRange 1 $totalNodes -dof 1 2 accel

# Record stresses for mesh elements using eleRange (recommended for Q4 elements)
recorder Element -eleRange 1 $totalElements -time -file RESULTS/eq_element_stress.out stress

# Create nodes file manually instead of using print command
set nodesFile [open "RESULTS/eq_nodes.txt" w]
puts $nodesFile "Node Information:"
for {set nodeID 1} {$nodeID <= $totalNodes} {incr nodeID} {
    puts $nodesFile "Node $nodeID: [nodeCoord $nodeID 1] [nodeCoord $nodeID 2]"
}
close $nodesFile

# Print only mesh elements (not boundary dashpots) - create custom file
set meshElementsFile [open "RESULTS/eq_elements.txt" w]
puts $meshElementsFile "Element Information:"
for {set j 0} {$j < $ny} {incr j} {
    for {set i 0} {$i < $nx} {incr i} {
        set ele_id [expr $j * $nx + $i + 1]
        set n1 [expr $i + $j*($nx+1) + 1]
        set n2 [expr $i + 1 + $j*($nx+1) + 1]
        set n3 [expr $i + 1 + ($j+1)*($nx+1) + 1]
        set n4 [expr $i + ($j+1)*($nx+1) + 1]
        puts $meshElementsFile "Element $ele_id: nodes $n1 $n2 $n3 $n4"
    }
}
close $meshElementsFile

# ----------------------
# Perform Analysis
# ----------------------
set dt_analysis 0.01
set t_final_analysis $tEnd
set num_steps [expr {int($t_final_analysis / $dt_analysis)}]

# Create model info file manually (moved here after dt_analysis is defined)
set modelInfoFile [open "RESULTS/eq_model_info.txt" w]
puts $modelInfoFile "Earthquake Model Information:"
puts $modelInfoFile "================================"
puts $modelInfoFile "Mesh Parameters:"
puts $modelInfoFile "  Width: $width m"
puts $modelInfoFile "  Height: $height m"
puts $modelInfoFile "  Elements: $nx x $ny = $totalElements"
puts $modelInfoFile "  Nodes: ($nx+1) x ($ny+1) = $totalNodes"
puts $modelInfoFile ""
puts $modelInfoFile "Material Properties:"
puts $modelInfoFile "  Young's Modulus (E): $E_soil Pa"
puts $modelInfoFile "  Poisson's Ratio (ν): $nu_soil"
puts $modelInfoFile "  Density (ρ): $rho_soil kg/m³"
puts $modelInfoFile "  Shear Wave Velocity (Vs): [format %.2f $Vs_soil] m/s"
puts $modelInfoFile "  P-Wave Velocity (Vp): [format %.2f $Vp_soil] m/s"
puts $modelInfoFile ""
puts $modelInfoFile "Analysis Parameters:"
puts $modelInfoFile "  Time Duration: $tEnd s"
puts $modelInfoFile "  Time Step: $dt_analysis s"
puts $modelInfoFile "  Earthquake Period: $period s"
puts $modelInfoFile "  Acceleration Factor: $accel_factor"
puts $modelInfoFile "  Rayleigh Damping (ζ): $zeta"
close $modelInfoFile

puts "Starting dynamic analysis (Q4 elements, $nx x $ny mesh)..."
analyze $num_steps $dt_analysis
puts "Dynamic analysis complete."