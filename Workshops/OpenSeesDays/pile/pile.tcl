###########################################################
#                                                         #
# Static pushover of a single pile, modeled as a beam on  #
#  a nonlinear Winkler foundation.  Lateral soil response #
#  is described by p-y springs.  Vertical soil response   #
#  described by t-z and q-z springs.                      #
#                                                         #
#   Created by:  Chris McGann                             #
#                HyungSuk Shin                            #
#                Pedro Arduino                            #
#                Peter Mackenzie-Helnwein                 #
#              --University of Washington--               #
#                                                         #
# ---> Basic units are kN and meters                      #
#                                                         #
###########################################################

#----------------------------------------------------------
#  pile geometry and mesh
#----------------------------------------------------------

set L1 1.0; # length of pile head (above ground surface) (m)
set L2 20.0; # length of embedded pile (below ground surface) (m)
set diameter 1.0; # pile diameter

set nElePile 84; # number of pile elements
set eleSize [expr ($L1+$L2)/$nElePile]; # pile element length 

# number of total pile nodes
set nNodePile [expr 1 + $nElePile]

#----------------------------------------------------------
#  create spring nodes
#----------------------------------------------------------

# spring nodes created with 3 dim, 3 dof
model BasicBuilder -ndm 3 -ndf 3

# counter to determine number of embedded nodes
set count 0

# create spring nodes
for {set i 1} {$i <= $nNodePile} {incr i} {
    set zCoord [expr $eleSize*($i - 1)]
  # only create spring nodes over embedded length of pile
    if {$zCoord <= $L2} {
        node $i             0.00   0.00   $zCoord
	node [expr $i+100]  0.00   0.00   $zCoord
        set count [expr $count+1]
    }
}

# number of embedded nodes
set nNodeEmbed $count

# spring node fixities
for {set i 1} {$i <= $nNodeEmbed} {incr i} {
    fix $i             1   1   1
    fix [expr $i+100]  0   1   1
}

#----------------------------------------------------------
#  soil properties
#----------------------------------------------------------

set gamma  17.0; # soil unit weight (kN/m^3)
set phi    36.0; # soil internal friction angle (degrees)
set Gsoil  150000; # soil shear modulus at pile tip (kPa)

set puSwitch 1 ; # select pult definition method for p-y curves 1: API 2 Brich
set kSwitch 1; # variation in coefficient of subgrade reaction with depth for p-y curves: 1 API, 2 modAPI
set gwtSwitch 1; # effect of ground water on subgrade reaction modulus for p-y curves: 1 above GWT, 2 below

#----------------------------------------------------------
#  create spring material objects
#----------------------------------------------------------

# specify necessary procedures
source get_pyParam.tcl
source get_tzParam.tcl
source get_qzParam.tcl

# p-y spring material
for {set i 1} {$i <= $nNodeEmbed} {incr i} {
  # depth of current py node
    set pyDepth [expr $L2 - $eleSize*($i - 1)]
  # procedure to define pult and y50
    set pyParam [get_pyParam $pyDepth $gamma $phi $diameter $eleSize $puSwitch $kSwitch $gwtSwitch]
    set pult [lindex $pyParam 0]
    set y50  [lindex $pyParam 1]
    uniaxialMaterial PySimple1 $i  2  $pult  $y50  0.0
}

# t-z spring material
for {set i 2} {$i <= $nNodeEmbed} {incr i} {
  # depth of current tz node
	set pyDepth [expr $eleSize*($i - 1)]
  # vertical effective stress at current depth
    set sigV    [expr $gamma*$pyDepth]
  # procedure to define tult and z50
    set tzParam  [get_tzParam $phi $diameter $sigV $eleSize]
    set tult [lindex $tzParam 0]
    set z50  [lindex $tzParam 1]
	uniaxialMaterial TzSimple1 [expr $i+100]  2  $tult  $z50  0.0
}

# q-z spring material
  # vertical effective stress at pile tip, no water table (depth is embedded pile length)
  set sigVq   [expr $gamma*$L2] 
  # procedure to define qult and z50
    set qzParam [get_qzParam $phi $diameter $sigVq $Gsoil]
    set qult [lindex $qzParam 0]
    set z50q [lindex $qzParam 1]
  
    uniaxialMaterial QzSimple1 101  2  $qult  $z50q  0.0  0.0 

puts "Finished creating all p-y, t-z, and z-z spring material objects..."

#----------------------------------------------------------
#  create zero-length elements for springs
#----------------------------------------------------------

# element at the pile tip (has q-z spring)
element zeroLength 1001  1  101  -mat 1  101  -dir 1  3

# remaining elements
for {set i 2} {$i <= $nNodeEmbed} {incr i} {
    element zeroLength [expr $i+1000]  $i  [expr $i+100]  -mat $i [expr $i+100]  -dir 1  3
}
puts "Finished creating all zero-Length elements for springs..."

#----------------------------------------------------------
#  create pile nodes
#----------------------------------------------------------

# pile nodes created with 3 dimensions, 6 degrees of freedom
model BasicBuilder -ndm 3 -ndf 6

# create pile nodes
for {set i 1} {$i <= $nNodePile} {incr i} {
   # z-coordinates of nodes depend on element length
    set zCoord [expr $eleSize*($i - 1)]

    node [expr $i+200] 0.00   0.00   $zCoord
}
puts "Finished creating all pile nodes..."

# create coordinate-transformation object
geomTransf Linear     1         0.0   -1.0  0.0

# create fixity at pile head (location of loading)
fix [expr 200+$nNodePile]     0    1    0    1    0    1

# create fixities for remaining pile nodes
for {set i 202} {$i < [expr 200+$nNodePile]} {incr i} {
    fix $i     0    1    0    1    0    1
}
puts "Finished creating all pile node fixities..."

#----------------------------------------------------------
#  define equal dof between pile and spring nodes
#----------------------------------------------------------

for {set i 1} {$i <= $nNodeEmbed} {incr i} {
   
    equalDOF [expr $i+200] [expr $i+100]  1  3
}
puts "Finished creating all equal degrees of freedom..."

#----------------------------------------------------------
#  pile section
#----------------------------------------------------------

# elastic pile section
source elasticPileSection.tcl

#----------------------------------------------------------
#  create pile elements
#----------------------------------------------------------

for {set i 201} {$i <= [expr 200+$nElePile]} {incr i} {
    element dispBeamColumn  $i  $i  [expr $i+1]  $secTag3D   3   1 
}
puts "Finished creating all pile elements..."

#----------------------------------------------------------
#  create recorders
#----------------------------------------------------------

# record information at specified increments
set timeStep 0.5

# record displacements at pile nodes
recorder Node -file pileDisp.out  -time  -nodeRange 201 [expr 200+$nNodePile]  -dof 1 2 3  -dT $timeStep  disp 
# record reaction force in the p-y springs
recorder Node -file reaction.out  -time  -nodeRange 1 $nNodePile  -dof 1  -dT $timeStep  reaction
# record element forces in pile elements
recorder Element -file pileForce.out  -time  -eleRange 201 [expr 200+$nElePile]  -dT $timeStep  globalForce
puts "Finished creating all recorders..."

# real time display recorder for visualization during analysis
recorder display "OpenSees Real Time" 10 10 700 700 -wipe
prp           0 100 0
vup           0 0 1
vpn           0 -1 0	
display       1 3 10

#----------------------------------------------------------
#  create the loading
#----------------------------------------------------------

setTime 10.0

# apply point load at the uppermost pile node in the x-direction
pattern Plain 10 {Series -time {0 10 20 10000} -values {0 0 1 1} -factor 1} {
	load [expr 200+$nNodePile]  3500  0.0  0.0  0.0  0.0  0.0
}
puts "Finished creating loading object..."

#----------------------------------------------------------
#  create the analysis
#----------------------------------------------------------

integrator LoadControl  0.05
numberer RCM
system SparseGeneral
constraints Transformation
test NormDispIncr 1e-5      20      1
algorithm Newton
analysis Static

set startT [clock seconds]
puts "Starting Load Application..."
analyze          201

set endT [clock seconds]
puts "Load Application finished..."
puts "Loading Analysis execution time: [expr $endT-$startT] seconds."

wipe
