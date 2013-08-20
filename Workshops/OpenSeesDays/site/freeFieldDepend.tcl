###########################################################
#                                                         #
# Site response analysis of a soil deposit with a         # 
# parabolically varying shear wave velocity profile on an #
# elastic half-space.  The finite rigidity of the         #
# underlying medium is considered through the use of a    #
# viscous damper at the base of the soil column.          #
#                                                         #
#   Created by:  Chris McGann                             #
#                HyungSuk Shin                            #
#                Pedro Arduino                            #
#                Peter Mackenzie-Helnwein                 #
#              --University of Washington--               #
#                                                         #
# ---> Basic units are kN and m   (unless specified)      #
#                                                         #
###########################################################

wipe

#-----------------------------------------------------------------------------------------------------------
#  1. DEFINE ANALYSIS PARAMETERS
#-----------------------------------------------------------------------------------------------------------

#---SOIL GEOMETRY
# number of soil layers
set numLayers       30
# layer thicknesses (m)
set layerThick(30)  1.0
set layerThick(29)  1.0
set layerThick(28)  1.0
set layerThick(27)  1.0
set layerThick(26)  1.0
set layerThick(25)  1.0
set layerThick(24)  1.0
set layerThick(23)  1.0
set layerThick(22)  1.0
set layerThick(21)  1.0
set layerThick(20)  1.0
set layerThick(19)  1.0
set layerThick(18)  1.0
set layerThick(17)  1.0
set layerThick(16)  1.0
set layerThick(15)  1.0
set layerThick(14)  1.0
set layerThick(13)  1.0
set layerThick(12)  1.0
set layerThick(11)  1.0
set layerThick(10)  1.0
set layerThick(9)   1.0
set layerThick(8)   1.0
set layerThick(7)   1.0
set layerThick(6)   1.0
set layerThick(5)   1.0
set layerThick(4)   1.0
set layerThick(3)   1.0
set layerThick(2)   1.0
set layerThick(1)   1.0

#---MATERIAL PROPERTIES
# soil mass density (Mg/m^3)
set rho             2.202
# soil shear wave velocity for each layer(m/s)
set Vs(30)          170.9
set Vs(29)          224.9
set Vs(28)          255.6
set Vs(27)          278.0
set Vs(26)          296.0
set Vs(25)          311.3
set Vs(24)          324.5
set Vs(23)          336.4
set Vs(22)          347.0
set Vs(21)          356.8
set Vs(20)          365.9
set Vs(19)          374.3
set Vs(18)          382.2
set Vs(17)          389.6
set Vs(16)          396.6
set Vs(15)          403.3
set Vs(14)          409.6
set Vs(13)          415.7
set Vs(12)          421.5
set Vs(11)          427.1
set Vs(10)          432.5
set Vs(9)           437.7
set Vs(8)           442.7
set Vs(7)           447.5
set Vs(6)           452.2
set Vs(5)           456.7
set Vs(4)           461.2
set Vs(3)           465.4
set Vs(2)           469.6
set Vs(1)           473.7
# soil shear modulus for each layer (kPa)
for {set k 1} {$k <= $numLayers} {incr k 1} {
    set G($k)       [expr $rho*$Vs($k)*$Vs($k)]
}
# poisson's ratio of soil
set nu              0.0
# soil elastic modulus for each layer (kPa)
for {set k 1} {$k <= $numLayers} {incr k 1} {
    set E($k)       [expr 2*$G($k)*(1+$nu)]
}
# soil bulk modulus for each layer (kPa)
for {set k 1} {$k <= $numLayers} {incr k 1} {
    set bulk($k)    [expr $E($k)/(3*(1-2*$nu))]
}
# soil friction angle
set phi             35.0
# peak shear strain
set gammaPeak       0.1
# reference pressure
set refPress        80.0
# pressure dependency coefficient
set pressCoeff      0.0
# phase transformation angle
set phaseAng        27.0
# contraction
set contract        0.06
# dilation coefficients
set dilate1         0.5
set dilate2         2.5
# liquefaction coefficients
set liq1            0.0 
set liq2            0.0 
set liq3            0.0
# bedrock shear wave velocity (m/s)
set rockVS           762.0
# bedrock mass density (Mg/m^3)
set rockDen          2.396

#---GROUND MOTION PARAMETERS
# time step in ground motion record
set motionDT        0.005
# number of steps in ground motion record
set motionSteps     7990

#---WAVELENGTH PARAMETERS
# highest frequency desired to be well resolved (Hz)
set fMax     100.0
# shear wave velocity desired to be well resolved
set vsMin    $Vs(30)
# wavelength of highest resolved frequency
set wave     [expr $vsMin/$fMax]
# number of elements per one wavelength
set nEle     8

#---RAYLEIGH DAMPING PARAMETERS
set pi      3.141592654
# damping ratio
set damp    0.02
# lower frequency
set omega1  [expr 2*$pi*0.2]
# upper frequency
set omega2  [expr 2*$pi*20]
# damping coefficients
set a0      [expr 2*$damp*$omega1*$omega2/($omega1 + $omega2)]
set a1      [expr 2*$damp/($omega1 + $omega2)]
puts "damping coefficients: a_0 = $a0;  a_1 = $a1"

#---ANALYSIS PARAMETERS
# Newmark parameters
set gamma           0.5
set beta            0.25

#-----------------------------------------------------------------------------------------------------------
#  2. DEFINE MESH GEOMETRY
#-----------------------------------------------------------------------------------------------------------

# trial vertical element size
set hTrial   [expr $wave/$nEle]

set numTotalEle 0
# loop over layers to determine appropriate element size
for {set k 1} {$k <= $numLayers} {incr k 1} {

  # trial number of elements 
    set nTrial   [expr $layerThick($k)/$hTrial]

  # round up if not an integer number of elements
    if { $nTrial > [expr floor($layerThick($k)/$hTrial)] } {
        set numEleY($k)  [expr int(floor($layerThick($k)/$hTrial)+1)]
    } else {
        set numEleY($k)  [expr int($nTrial)]
    }
    puts "number of vertical elements in layer $k: $numEleY($k)"

  # counter for total number of elements
    set numTotalEle [expr $numTotalEle + $numEleY($k)]

  # final vertical size of elements (m) 
    set sizeEleY($k)  [expr {$layerThick($k)/$numEleY($k)}]
    puts "vertical size of elements in layer $k: $sizeEleY($k)"
}
puts "total number of vertical elements: $numTotalEle"

# define horizontal size of elements as smallest vertical element size (m)
set sizeEleX  $sizeEleY(1)
for {set k 2} {$k <= $numLayers} {incr k 1} {
    if { $sizeEleY($k) < $sizeEleY([expr $k-1]) } {
        set sizeEleX  $sizeEleY($k)
    }
}
puts "horizontal size of elements: $sizeEleX"

# number of nodes in vertical direction in each layer
set numTotalNode 0
for {set k 1} {$k < $numLayers} {incr k 1} {
    set numNodeY($k)  [expr 2*$numEleY($k)]
    puts "number of nodes in layer $k: $numNodeY($k)"
    set numTotalNode  [expr $numTotalNode + $numNodeY($k)]
}
set numNodeY($numLayers) [expr 2*($numEleY($numLayers)+1)]
puts "number of nodes in layer $numLayers: $numNodeY($numLayers)"
set numTotalNode      [expr $numTotalNode + $numNodeY($numLayers)]
puts "total number of nodes: $numTotalNode"

#-----------------------------------------------------------------------------------------------------------
#  3. DEFINE NODES FOR SOIL ELEMENTS
#-----------------------------------------------------------------------------------------------------------

# soil nodes are created in 2 dimensions, with 3 dof (2 translational, 1 porePressure)
model BasicBuilder -ndm 2 -ndf 2

set yCoord     0.0
set count      0
# loop over layers
for {set k 1} {$k <= $numLayers} {incr k 1} {
  # loop over nodes
    for {set j 1} {$j <= $numNodeY($k)} {incr j 2} {
        node    [expr $j+$count]     0.0         $yCoord
        node    [expr $j+$count+1]   $sizeEleX   $yCoord

        set yCoord  [expr {$yCoord + $sizeEleY($k)}]
    }
    set count [expr $count+$numNodeY($k)]
}
puts "Finished creating all soil nodes..."

#-----------------------------------------------------------------------------------------------------------
#  4. DEFINE DASHPOT NODES
#-----------------------------------------------------------------------------------------------------------

node 2000 0.0 0.0
node 2001 0.0 0.0

puts "Finished creating dashpot nodes..."

#-----------------------------------------------------------------------------------------------------------
#  5. DEFINE BOUNDARY CONDITIONS AND EQUAL DOF
#-----------------------------------------------------------------------------------------------------------

# define fixity of base nodes
fix 1 0 1
fix 2 0 1

# define fixity of dashpot nodes
fix  2000  1 1
fix  2001  0 1

# define equal DOF for simple shear deformation of soil elements
for {set k 3} {$k <= $numTotalNode} {incr k 2} {
    equalDOF  $k  [expr $k+1]  1 2
}

# define equal DOF for dashpot and base soil nodes
equalDOF 1    2 1 
equalDOF 1 2001 1

puts "Finished creating all boundary conditions and equalDOF..."

#-----------------------------------------------------------------------------------------------------------
#  6. DEFINE SOIL MATERIALS
#-----------------------------------------------------------------------------------------------------------

# loop over layers to define materials
for {set k 1} {$k <= $numLayers} {incr k 1} {
nDMaterial PressureDependMultiYield $k 2 $rho $G($k) $bulk($k) $phi $gammaPeak $refPress $pressCoeff \
                                    $phaseAng $contract $dilate1 $dilate2 $liq1 $liq2 $liq3 16
}

puts "Finished creating all soil materials..."

#-----------------------------------------------------------------------------------------------------------
#  7. DEFINE SOIL ELEMENTS
#-----------------------------------------------------------------------------------------------------------

set wgtX  0.0
set wgtY  [expr -9.81*$rho]

set count 0
# loop over layers 
for {set k 1} {$k <= $numLayers} {incr k 1} {
  # loop over elements
    for {set j 1} {$j <= $numEleY($k)} {incr j 1} {

        set nI  [expr 2*($j+$count) - 1] 
        set nJ  [expr $nI + 1]
        set nK  [expr $nI + 3]
        set nL  [expr $nI + 2]

        element quad [expr $j+$count] $nI $nJ $nK $nL 1.0 "PlaneStrain" $k 0.0 0.0 $wgtX $wgtY

    }
    set count [expr $count + $numEleY($k)]
}

puts "Finished creating all soil elements..."

#-----------------------------------------------------------------------------------------------------------
#  8. DEFINE MATERIAL AND ELEMENTS FOR VISCOUS DAMPERS
#-----------------------------------------------------------------------------------------------------------

# dashpot coefficient 
set mC     [expr $sizeEleX*$rockDen*$rockVS]

# material
uniaxialMaterial Viscous 4000 $mC 1

# elements
element zeroLength 5000 2000 2001 -mat 4000 -dir 1

puts "Finished creating dashpot material and element..."

#-----------------------------------------------------------------------------------------------------------
#  8. CREATE GRAVITY RECORDERS
#-----------------------------------------------------------------------------------------------------------

# record nodal displacments, velocities, and accelerations at each time step
recorder Node -file Gdisplacement.out -time  -nodeRange 1 $numTotalNode -dof 1 2  disp
recorder Node -file Gvelocity.out     -time  -nodeRange 1 $numTotalNode -dof 1 2  vel
recorder Node -file Gacceleration.out -time  -nodeRange 1 $numTotalNode -dof 1 2  accel

# record stress and strain at each gauss point in the soil elements
recorder Element -file Gstress1.out   -time  -eleRange  1   $numTotalEle   material 1 stress
recorder Element -file Gstress2.out   -time  -eleRange  1   $numTotalEle   material 2 stress
recorder Element -file Gstress3.out   -time  -eleRange  1   $numTotalEle   material 3 stress
recorder Element -file Gstress4.out   -time  -eleRange  1   $numTotalEle   material 4 stress

recorder Element -file Gstrain1.out   -time  -eleRange  1   $numTotalEle   material 1 strain
recorder Element -file Gstrain2.out   -time  -eleRange  1   $numTotalEle   material 2 strain
recorder Element -file Gstrain3.out   -time  -eleRange  1   $numTotalEle   material 3 strain
recorder Element -file Gstrain4.out   -time  -eleRange  1   $numTotalEle   material 4 strain

puts "Finished creating gravity recorders..."

#-----------------------------------------------------------------------------------------------------------
#  9. APPLY GRAVITY LOADING
#-----------------------------------------------------------------------------------------------------------

constraints Transformation
test        NormDispIncr 1e-5 30 1
algorithm   Newton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta 
analysis    Transient

analyze     10 5.0e2

puts "Finished with elastic gravity analysis..."

# update materials to consider plastic behavior
for {set k 1} {$k <= $numLayers} {incr k} {
    updateMaterialStage -material $k -stage 1
}

# plastic gravity loading
analyze     40 5.0e2

puts "Finished with plastic gravity analysis..."

#-----------------------------------------------------------------------------------------------------------
#  10. CREATE POST-GRAVITY RECORDERS
#-----------------------------------------------------------------------------------------------------------

# reset time and analysis
setTime 0.0
wipeAnalysis

# record nodal displacments, velocities, and accelerations at each time step
recorder Node -file displacement.out -time -dT $motionDT -nodeRange 1 $numTotalNode -dof 1 2  disp
recorder Node -file velocity.out     -time -dT $motionDT -nodeRange 1 $numTotalNode -dof 1 2  vel
recorder Node -file acceleration.out -time -dT $motionDT -nodeRange 1 $numTotalNode -dof 1 2  accel

# record stress and strain at each gauss point in the soil elements
recorder Element -file stress1.out   -time -dT $motionDT -eleRange  1   $numTotalEle   material 1 stress
recorder Element -file stress2.out   -time -dT $motionDT -eleRange  1   $numTotalEle   material 2 stress
recorder Element -file stress3.out   -time -dT $motionDT -eleRange  1   $numTotalEle   material 3 stress
recorder Element -file stress4.out   -time -dT $motionDT -eleRange  1   $numTotalEle   material 4 stress

recorder Element -file strain1.out   -time -dT $motionDT -eleRange  1   $numTotalEle   material 1 strain
recorder Element -file strain2.out   -time -dT $motionDT -eleRange  1   $numTotalEle   material 2 strain
recorder Element -file strain3.out   -time -dT $motionDT -eleRange  1   $numTotalEle   material 3 strain
recorder Element -file strain4.out   -time -dT $motionDT -eleRange  1   $numTotalEle   material 4 strain

puts "Finished creating all recorders..."

#-----------------------------------------------------------------------------------------------------------
#  11. DYNAMIC ANALYSIS  
#-----------------------------------------------------------------------------------------------------------

# define constant factor for applied velocity
set cFactor [expr $sizeEleX*$rockDen*$rockVS]

# define velocity time history file
set velocityFile velocityHistory.out

# timeseries object for force history
set mSeries "Path -dt $motionDT -filePath $velocityFile -factor $cFactor"

# loading object
pattern Plain 10 $mSeries {
load 1 1.0 0.0 
}
puts "Dynamic loading created..."

#---DETERMINE STABLE ANALYSIS TIME STEP USING CFL CONDITION
# maximum shear wave velocity (m/s)
set vsMax     $Vs(1)
# duration of ground motion (s)
set duration  [expr $motionDT*$motionSteps]
# trial analysis time step
set kTrial    [expr $sizeEleX/(pow($vsMax,0.5))]
# define time step and number of steps for analysis
if { $motionDT <= $kTrial } {
    set nSteps  $motionSteps
    set dT      $motionDT
} else {
    set nSteps  [expr floor($duration/$kTrial)+1]
    set dT      [expr $duration/$nSteps] 
}
puts "number of steps in analysis: $nSteps"
puts "analysis time step: $dT"

# analysis objects
constraints Transformation
test        NormDispIncr 1e-3 35 1
algorithm   Newton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta 
rayleigh    $a0 $a1 0.0 0.0
analysis    Transient

analyze     $nSteps  $dT

puts "Finished with dynamic analysis..."

wipe
