#############################################################################
#-------Department of Civil, Structural and Environmental Engineering-------#
#---------------------------University at Buffalo---------------------------#
# Modeling of Triple FP isolator                                            #
# Written By: Hyun-myung Kim (hkim59@buffalo.edu)                           #
# Date: May, 2023                                                           #
#############################################################################

# Units: N, m, sec
# Remove existing model
wipe

# EXAMPLE 3 (Kim and Constantinou, 2023 https://doi.org/10.1002/eqe.3797)
#----------------------------------------------------------------------------
# User Defined Parameters
#----------------------------------------------------------------------------

# TFP Geomoetry of Configuration A
set L1 0.3937;                    # Effective radii (m)
set L2 3.7465;
set L3 3.7465;
set d1 0.0716;                    # Actual displacement capacity (m)
set d2 0.5043;
set d3 0.5043;
set b1 [expr 0.508];              # Diameter of the rigid slider and the two inner slide plate (m)
set b2 [expr 0.711];
set b3 [expr 0.711];
set r1 [expr $b1/2];              # Radius of of the rigid slider and the two inner slide plate (m)
set r2 [expr $b2/2];
set r3 [expr $b3/2];

set uy 0.001;                     # Yield displacement (m)
set kvc 8000000000.;              # Vertical compression stiffness (N/m)
set kvt 1.;                       # Vertical tension stiffness (N/m)
set minFv 0.1;                    # Minimum compression force in the bearing (N)

set g     9.81;                   # Gravity acceleration (m/s^2)
set P     13345e+03;              # Load on top of TFP
set Mass [expr $P/$g];            # Mass on top of TFP
set tol 1.e-5;                    # Relative tolerance for checking convergence

# Heat parameters
set Diffu 0.444e-5;               # Thermal diffusivity (m^2/sec)
set Conduct 18;                   # Thermal conductivity (W/m*Celsius)
set Temperature0 20;              # Initial temperature (Celsius)

# Friction coefficients (reference)
set mu1 0.01;
set mu2 0.04;
set mu3 0.08;

# Reference Pressure
set Pref1 [expr $P/($r1*$r1*3.141592)];
set Pref2 [expr $P/($r2*$r2*3.141592)];
set Pref3 [expr $P/($r3*$r3*3.141592)];

#----------------------------------------------------------------------------
# Start of model generation
#----------------------------------------------------------------------------

#Create Model Builder
model basic -ndm 3 -ndf 6

# Create nodes
node 1 0 0 0; # End i
node 2 0 0 0; # End j

# Define single point constraints
fix 1     1 1 1 1 1 1;

# Define friction models
set tagTemp 1;
set tagVel 0;
set tagPres 0;
set velRate 100;
set kTmodel 1;                     # kT = 1/2 at 200 degree celsius

#----------------------------------------------------------------------------
# Bring material models and define element
#----------------------------------------------------------------------------

# Creating material for compression and rotation behaviors
uniaxialMaterial Elastic 1 $kvc;
uniaxialMaterial Elastic 2 10.;

set tagT 1;

# Define TripleFrictionPendulumX element
# element TripleFrictionPendulumX $eleTag $iNode $jNode $tagT $vertMatTag $rotZMatTag $rotXMatTag $rotYMatTag $tagPres $tagTemp $tagVel $mu1 $mu2 $mu3 $L1 $L2 $L3 $d1 $d2 $d3 $b1 $b2 $b3 $W $uy $kvt $minFv $tol $Pref1 $Pref2 $Pref3 $Diffu $Conduct $Temperature0 $velRate $kTmodel $unit
element TripleFrictionPendulumX 1 1 2  $tagT  1 2 2 2 $tagPres $tagTemp $tagVel $mu1 $mu2 $mu3 $L1 $L2 $L3 $d1 $d2 $d3 $b1 $b2 $b3 $P $uy $kvt $minFv $tol $Pref1 $Pref2 $Pref3 $Diffu $Conduct $Temperature0 $velRate $kTmodel 1;

#----------------------------------------------------------------------------
# Apply gravity load
#----------------------------------------------------------------------------

#Create a plain load pattern with linear timeseries
pattern Plain 1 "Linear" {

        load 2 0. 0. -[expr $P] 0.0 0.0 0.0
}

#----------------------------------------------------------------------------
# Start of analysis generation (Gravity)
#----------------------------------------------------------------------------

system BandSPD
constraints Transformation
numberer RCM
test NormDispIncr 1.0e-15 10 3
algorithm Newton
integrator LoadControl 1
analysis Static

#----------------------------------------------------------------------------
# Analysis (Gravity)
#----------------------------------------------------------------------------

analyze 1
puts "Gravity analysis completed SUCCESSFULLY";

#----------------------------------------------------------------------------
# Start of analysis generation
# (Sinusoidal; Two cycles of 5s period and 508mm amplitude)
#----------------------------------------------------------------------------

loadConst -time 0.0

# analysis time step
set dt [expr 0.008]

# excitation time step
set dt1 [expr 0.001]

timeSeries Trig 11 $dt 10 5 -factor 0.508 -shift 0

pattern MultiSupport 2 {
groundMotion 1 Plain -disp 11
# Node, direction, GMtag
imposedMotion 2 2 1
}

#----------------------------------------------------------------------------
# Start of recorder generation (Sinusoidal)
#----------------------------------------------------------------------------

# Set up recorder
set OutDir                EXAMPLE3;                       # Output folder
set OutFile1      TEMPERATURE.txt;
set OutFile2      DISP.txt;
set OutFile3      FORCE.txt;
set OutFile4      COMPDISP.txt;

file mkdir $OutDir;
recorder Element -file $OutDir/$OutFile1 -time -ele 1 Parameters;
recorder Node -file $OutDir/$OutFile2 -time -nodes 2 -dof 1 2 3 disp;
recorder Element -file $OutDir/$OutFile3 -time -ele 1 basicForce;
recorder Element -file $OutDir/$OutFile4 -time -ele 1 compDisplacement;

#----------------------------------------------------------------------------
# Analysis (Sinusoidal)
#----------------------------------------------------------------------------

system SparseGeneral
constraints Transformation
test NormDispIncr 1.0e-5 20 0
algorithm Newton
numberer Plain
integrator Newmark 0.5 0.25
analysis Transient

# set some variables
set tFinal [expr 10]
set tCurrent [getTime]
set ok 0

# Perform the transient analysis
while {$ok == 0 && $tCurrent < $tFinal} {

    set ok [analyze 1 $dt]

# if the analysis fails try initial tangent iteration
  if {$ok != 0} {
  puts "regular newton failed .. lets try an initail stiffness for this step"
  test NormDispIncr 1.0e-12  100 0
  algorithm ModifiedNewton -initial
  set ok [analyze 1 $dt]
  if {$ok == 0} {puts "that worked .. back to regular newton"}
  test NormDispIncr 1.0e-12  10
  algorithm Newton
  }

    set tCurrent [getTime]
}

# Print a message to indicate if analysis succesfull or not
if {$ok == 0} {
   puts "Transient analysis completed SUCCESSFULLY";
} else {
   puts "Transient analysis completed FAILED";
}