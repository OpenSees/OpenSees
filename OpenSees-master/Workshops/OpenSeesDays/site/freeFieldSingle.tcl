###########################################################
#                                                         #
# Site response analysis of a soil deposit on an elastic  #
# half-space.  The finite rigidity of the underlying      #
# medium is considered through the use of a viscous       #
# damper at the base of the soil column.                  #
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

#-------------------------------------------------------------------------------------------
#  1. DEFINE ANALYSIS PARAMETERS
#-------------------------------------------------------------------------------------------

#---SOIL GEOMETRY
# thicknesses of soil profile (m)
set soilThick  40.0

#---MATERIAL PROPERTIES
# soil mass density (Mg/m^3)
set rho         1.7
# soil shear wave velocity (m/s)
set Vs          250.0
# soil shear modulus (kPa)
set G           [expr $rho*$Vs*$Vs]
# poisson's ratio of soil
set nu          0.0
# soil elastic modulus (kPa)
set E           [expr 2*$G*(1+$nu)]
# soil bulk modulus (kPa)
set bulk        [expr $E/(3*(1-2*$nu))]
# soil cohesion (kPa)
set cohesion    95.0
# peak shear strain
set gammaPeak   0.05
# soil friction angle
set phi         0.0
# reference pressure
set refPress    80.0
# pressure dependency coefficient
set pressCoeff  0.0
# bedrock shear wave velocity (m/s)
set rockVS      760
# bedrock mass density (Mg/m^3)
set rockDen     2.4

#---GROUND MOTION PARAMETERS
# time step in ground motion record
set motionDT        0.005
# number of steps in ground motion record
set motionSteps     7990

#---WAVELENGTH PARAMETERS
# highest frequency desired to be well resolved (Hz)
set fMax     100.0
# wavelength of highest resolved frequency
set wave     [expr $Vs/$fMax]
# number of elements per one wavelength
set nEle     10

#---RAYLEIGH DAMPING PARAMETERS
set pi      3.141592654
# damping ratio
set damp    0.02
# lower frequency
set omega1     [expr 2*$pi*0.2]
# upper frequency
set omega2     [expr 2*$pi*20]
# damping coefficients
set a0      [expr 2*$damp*$omega1*$omega2/($omega1 + $omega2)]
set a1      [expr 2*$damp/($omega1 + $omega2)]
puts "damping coefficients: a_0 = $a0;  a_1 = $a1"

#---ANALYSIS PARAMETERS
# Newmark parameters
set gamma           0.5
set beta            0.25

#-------------------------------------------------------------------------------------------
#  2. DEFINE MESH GEOMETRY
#-------------------------------------------------------------------------------------------

# trial vertical element size
set hTrial    [expr $wave/$nEle]

# trial number of elements 
set nTrial    [expr $soilThick/$hTrial]
# round up if not an integer number of elements
if { $nTrial > [expr floor($soilThick/$hTrial)] } {
    set numEleY  [expr int(floor($soilThick/$hTrial)+1)]
} else {
    set numEleY  [expr int($nTrial)]
}
puts "number of verticl elements: $numEleY"

# final vertical size of elements (m) 
set sizeEleY  [expr $soilThick/$numEleY]
puts "vertical size of elements: $sizeEleY"

# define horizontal size of elements (m)
set sizeEleX  $sizeEleY
puts "horizontal size of elements: $sizeEleX"

# number of nodes in vertical direction
set numNodeY  [expr 2*($numEleY + 1)]

#-------------------------------------------------------------------------------------------
#  3. DEFINE NODES FOR SOIL ELEMENTS
#-------------------------------------------------------------------------------------------

# soil nodes are created in 2 dimensions, with 2 translational dof
model BasicBuilder -ndm 2 -ndf 2

set yCoord     0.0
set count      0
# loop over nodes
for {set j 1} {$j <= $numNodeY} {incr j 2} {
    node    $j            0.0         [expr $yCoord + $count*$sizeEleY]
    node    [expr $j+1]   $sizeEleX   [expr $yCoord + $count*$sizeEleY]

    set count [expr $count+1]
}
puts "Finished creating all soil nodes..."

#-------------------------------------------------------------------------------------------
#  4. DEFINE DASHPOT NODES
#-------------------------------------------------------------------------------------------

node 2000 0.0 0.0
node 2001 0.0 0.0

puts "Finished creating dashpot nodes..."

#-------------------------------------------------------------------------------------------
#  5. DEFINE BOUNDARY CONDITIONS AND EQUAL DOF
#-------------------------------------------------------------------------------------------

# define fixity of base nodes
fix 1 0 1
fix 2 0 1

# define fixity of dashpot nodes
fix  2000  1 1
fix  2001  0 1

# define equal DOF for simple shear deformation of soil elements
for {set k 3} {$k <= $numNodeY} {incr k 2} {
    equalDOF  $k  [expr $k+1]  1 2
}

# define equal DOF for dashpot and base soil nodes
equalDOF 1    2 1 
equalDOF 1 2001 1

puts "Finished creating all boundary conditions and equalDOF..."

#-------------------------------------------------------------------------------------------
#  6. DEFINE SOIL MATERIALS
#-------------------------------------------------------------------------------------------

nDMaterial PressureIndependMultiYield 1 2 $rho $G $bulk $cohesion $gammaPeak \
                                          $phi $refPress $pressCoeff 25

puts "Finished creating all soil materials..."

#-------------------------------------------------------------------------------------------
#  7. DEFINE SOIL ELEMENTS
#-------------------------------------------------------------------------------------------

set wgtX  0.0
set wgtY  [expr -9.81*$rho]

# loop over elements
for {set j 1} {$j <= $numEleY} {incr j 1} {

        set nI  [expr 2*$j - 1] 
        set nJ  [expr $nI + 1]
        set nK  [expr $nI + 3]
        set nL  [expr $nI + 2]

        element quad $j $nI $nJ $nK $nL 1.0 "PlaneStrain" 1 0.0 0.0 $wgtX $wgtY
}

puts "Finished creating all soil elements..."

#-------------------------------------------------------------------------------------------
#  8. DEFINE MATERIAL AND ELEMENTS FOR VISCOUS DAMPERS
#-------------------------------------------------------------------------------------------

# dashpot coefficient 
set mC     [expr $sizeEleX*$rockDen*$rockVS]

# material
uniaxialMaterial Viscous 4000 $mC 1

# elements
element zeroLength 5000 2000 2001 -mat 4000 -dir 1

puts "Finished creating dashpot material and element..."

#-------------------------------------------------------------------------------------------
#  8. CREATE GRAVITY RECORDERS
#-------------------------------------------------------------------------------------------

# record nodal displacements, velocities, and accelerations at each time step
recorder Node -file Gdisplacement.out -time  -nodeRange 1 $numNodeY -dof 1 2  disp
recorder Node -file Gvelocity.out     -time  -nodeRange 1 $numNodeY -dof 1 2  vel
recorder Node -file Gacceleration.out -time  -nodeRange 1 $numNodeY -dof 1 2  accel

# record stress and strain at each gauss point in the soil elements
recorder Element -file Gstress1.out   -time  -eleRange  1   $numEleY   material 1 stress
recorder Element -file Gstress2.out   -time  -eleRange  1   $numEleY   material 2 stress
recorder Element -file Gstress3.out   -time  -eleRange  1   $numEleY   material 3 stress
recorder Element -file Gstress4.out   -time  -eleRange  1   $numEleY   material 4 stress

recorder Element -file Gstrain1.out   -time  -eleRange  1   $numEleY   material 1 strain
recorder Element -file Gstrain2.out   -time  -eleRange  1   $numEleY   material 2 strain
recorder Element -file Gstrain3.out   -time  -eleRange  1   $numEleY   material 3 strain
recorder Element -file Gstrain4.out   -time  -eleRange  1   $numEleY   material 4 strain

puts "Finished creating gravity recorders..."

#-------------------------------------------------------------------------------------------
#  9. APPLY GRAVITY LOADING
#-------------------------------------------------------------------------------------------

constraints Transformation
test        NormDispIncr 1e-5 30 1
algorithm   Newton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta 
analysis    Transient

analyze     10 5.0e2

puts "Finished with elastic gravity analysis..."

# update material to consider elastoplastic behavior
updateMaterialStage -material 1 -stage 1

# plastic gravity loading
analyze     40 5.0e2

puts "Finished with plastic gravity analysis..."

#-------------------------------------------------------------------------------------------
#  10. CREATE POST-GRAVITY RECORDERS
#-------------------------------------------------------------------------------------------

# reset time and analysis
setTime 0.0
wipeAnalysis

# record nodal displacements, velocities, and accelerations at each time step
recorder Node -file displacement.out -time -dT $motionDT -nodeRange 1 $numNodeY -dof 1 2  disp
recorder Node -file velocity.out     -time -dT $motionDT -nodeRange 1 $numNodeY -dof 1 2  vel
recorder Node -file acceleration.out -time -dT $motionDT -nodeRange 1 $numNodeY -dof 1 2  accel

# record stress and strain at each gauss point in the soil elements
recorder Element -file stress1.out   -time -dT $motionDT -eleRange  1   $numEleY   material 1 stress
recorder Element -file stress2.out   -time -dT $motionDT -eleRange  1   $numEleY   material 2 stress
recorder Element -file stress3.out   -time -dT $motionDT -eleRange  1   $numEleY   material 3 stress
recorder Element -file stress4.out   -time -dT $motionDT -eleRange  1   $numEleY   material 4 stress

recorder Element -file strain1.out   -time -dT $motionDT -eleRange  1   $numEleY   material 1 strain
recorder Element -file strain2.out   -time -dT $motionDT -eleRange  1   $numEleY   material 2 strain
recorder Element -file strain3.out   -time -dT $motionDT -eleRange  1   $numEleY   material 3 strain
recorder Element -file strain4.out   -time -dT $motionDT -eleRange  1   $numEleY   material 4 strain

puts "Finished creating all recorders..."

#-------------------------------------------------------------------------------------------
#  11. DYNAMIC ANALYSIS  
#-------------------------------------------------------------------------------------------

# define constant factor for applied velocity
set cFactor [expr $sizeEleX*$rockDen*$rockVS]

# define velocity time history file
set velocityFile velocityHistory.out

# timeseries object for applied force history
set mSeries "Path -dt $motionDT -filePath $velocityFile -factor $cFactor"

# loading object
pattern Plain 10 $mSeries {
    load 1 1.0 0.0 
}
puts "Dynamic loading created..."

#---DETERMINE STABLE ANALYSIS TIME STEP USING CFL CONDITION
# duration of ground motion (s)
set duration  [expr $motionDT*$motionSteps]
# trial analysis time step 
set kTrial    [expr $sizeEleX/(pow($Vs,0.5))]
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
test        NormDispIncr 1e-3 15 1
algorithm   Newton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta 
rayleigh    $a0 $a1 0.0 0.0
analysis    Transient

analyze     $nSteps  $dT

puts "Finished with dynamic analysis..."

wipe
