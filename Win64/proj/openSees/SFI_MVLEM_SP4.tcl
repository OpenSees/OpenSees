# ----------------------------------------------------------------
# Example 1: Simulation of wall cyclic behavior using SFI_MVLEM
# Specimen: RW-A15-P10-S78 (Tran and Wallace, 2012)
# Created by: Kristijan Kolozvari (kkolozvari@fullerton.edu)
# Date: 7/2015
# ----------------------------------------------------------------

# --------------------------------------------------------
# Start of model generation (Units: kip, in, sec, ksi)
# --------------------------------------------------------
wipe 

# Increment
set Dincr 0.01; # displacement increment for displacement controlled analysis. 
# Shear resisting mechanism parameters
set nu 1.0; 		# friction coefficient
set alfadow 0.01; 	# dowel action stiffness parameter

# Set Up Directories
#set modelName "SP4"; 				# Model Name
#set modelName "$Dincr _ $alfadow _ $nu"; 				# Model Name
#set dataDir SFI_MVLEM_$modelName;		# Name of output folder
#file mkdir $dataDir;

# Create ModelBuilder for 2D element (with two-dimensions and 2 DOF/node)
model BasicBuilder -ndm 2 -ndf 3

# --------------------------------------------
# Set geometry, nodes, boundary conditions
# --------------------------------------------

# Wall Geometry 
set H 72; 	# Wall height
set t 6;	# Wall thickness

# Create nodes
# node nodeId xCrd yCrd
node 1 0 0;  
node 2 0 12;
node 3 0 24;  
node 4 0 40;  
node 5 0 56;  
node 6 0 $H;
 
# Boundary conditions
fix 1 1 1 1;	# Fixed condition at node 1
 
# Set Control Node and DOF
set IDctrlNode 6;
set IDctrlDOF 1;

# ------------------------------------------------------------------------
# Define uniaxial materials for 2D RC Panel Constitutive Model (FSAM)
# ------------------------------------------------------------------------

# STEEL ...........................................................
# uniaxialMaterial SteelMPF $mattag  $fyp  $fyn  $E0  $bp  $bn   $R0  $a1  $a2
# steel X
set fyX 58.4103;	# fy
set bx 0.01; 		# strain hardening

# steel Y web
set fyYw 58.4103; 	# fy
set byw 0.01; 		# strain hardening

# steel Y boundary
set fyYb 69.0; 		# fy
set byb 0.002; 		# strain hardening

# steel misc
set Esy 29000.0; 	# Young's modulus
set Esx $Esy; 		# Young's modulus
set R0 20.0; 		# initial value of curvature parameter
set A1 0.925;		# curvature degradation parameter
set A2 0.15;		# curvature degradation parameter

# Build steel materials
uniaxialMaterial    SteelMPF  1 $fyX $fyX $Esx $bx $bx $R0 $A1 $A2; 	# steel X
uniaxialMaterial    SteelMPF  2 $fyYw $fyYw $Esy $byw $byw $R0 $A1 $A2; # steel Y web
uniaxialMaterial    SteelMPF  3 $fyYb $fyYb $Esy $byb $byb $R0 $A1 $A2; # steel Y boundary

# CONCRETE ........................................................
# uniaxialMaterial ConcreteCM $mattag  $fpcc  $epcc  $Ec  $rc  $xcrn   $ft  $et  $rt   $xcrp <-GapClose $gap>

# unconfined
set fpc 8.09; 		# peak compressive stress
set ec0 -0.002371;	# strain at peak compressive stress
set ft 0.335798;	# peak tensile stress
set et 0.00008;		# strain at peak tensile stress
set Ec 5403.2172; 	# Young's modulus	
set xcrnu 1.022;	# cracking strain - compression
set xcrp 10000;		# cracking strain - tension
set ru 15;			# shape parameter - compression
set rt 1.2;			# shape parameter - tension

# confined
set fpcc 10.479723; # peak compressive stress
set ec0c -0.005873; # strain at peak compressive stress
set Ecc 5953.9187; 	# Young's modulus
set xcrnc 1.023;	# cracking strain - compression
set rc 12.072964;	# shape parameter - compression

# Build concrete materials
#uniaxialMaterial ConcreteCM 4 -$fpc   $ec0   $Ec  $ru  $xcrnu  $ft  $et  $rt  $xcrp -GapClose 0; # unconfined concrete
#uniaxialMaterial ConcreteCM 5 -$fpcc  $ec0c  $Ecc $rc  $xcrnc  $ft  $et  $rt  $xcrp -GapClose 0; # confined concrete

# uniaxialMaterial Concrete02 $matTag $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets
uniaxialMaterial Concrete02 4 -$fpc   $ec0   [expr -0.2*$fpc] -0.015 0.3 $ft  300; # unconfined concrete
uniaxialMaterial Concrete02 5 -$fpcc  $ec0c  [expr -0.2*$fpc] -0.10  0.3 $ft  300; # confined concrete

# ---------------------------------------
#  Define 2D RC Panel Material (FSAM)
# ---------------------------------------

# Reinforcing ratios
set rouXw 0.0074; 	# X web
set rouXb 0.0082; 	# X boundary
set rouYw 0.0074; 	# Y web
set rouYb 0.0587; 	# Y boundary

# nDMaterial FSAM $mattag   $rho  $sX  $sY  $conc  $rouX  $rouY  $nu  $alfadow
nDMaterial FSAM         6    0.0   1     2     4  $rouXw $rouYw  $nu  $alfadow; # Web (unconfined concrete)
nDMaterial FSAM         7    0.0   1     3     5  $rouXb $rouYb  $nu  $alfadow; # Boundary (confined concrete only)

# ------------------------------
#  Define SFI_MVLEM elements
# ------------------------------

# element SFI_MVLEM eleTag iNode jNode m c -thick fiberThick -width fiberWidth -mat matTags 
element SFI_MVLEM  1  1  2 5 0.4 -thick $t $t $t $t $t -width 9 10 10 10 9 -mat 7 6 6 6 7
element SFI_MVLEM  2  2  3 5 0.4 -thick $t $t $t $t $t -width 9 10 10 10 9 -mat 7 6 6 6 7
element SFI_MVLEM  3  3  4 5 0.4 -thick $t $t $t $t $t -width 9 10 10 10 9 -mat 7 6 6 6 7
element SFI_MVLEM  4  4  5 5 0.4 -thick $t $t $t $t $t -width 9 10 10 10 9 -mat 7 6 6 6 7
element SFI_MVLEM  5  5  6 5 0.4 -thick $t $t $t $t $t -width 9 10 10 10 9 -mat 7 6 6 6 7

# ------------------------------
# End of model generation
# ------------------------------

# Initialize
initialize

# ------------------------------
# Recorder generation
# ------------------------------

# Nodal recorders
recorder Node -file MVLEM_Dtop.out -time -node $IDctrlNode -dof 1 disp
#recorder Node -file $dataDir/MVLEM_DOFs.out -time -node 1 2 3 4 5 6 -dof 1 2 3 disp

# Element recorders
recorder Element -file MVLEM_Fgl.out -time -ele 1 2 3 4 5 globalForce
recorder Element -file MVLEM_Dsh.out -time -ele 1 2 3 4 5 ShearDef
recorder Element -file MVLEM_Curvature.out -time -ele 1 2 3 4 5 Curvature

# Single RC panel (macro-fiber) responses
#recorder Element -file $dataDir/MVLEM_panel_strain.out -time -ele 1 RCPanel 1 panel_strain
#recorder Element -file $dataDir/MVLEM_panel_stress.out -time -ele 1 RCPanel 1 panel_stress
#recorder Element -file $dataDir/MVLEM_panel_stress_concrete.out -time -ele 1 RCPanel 1 panel_stress_concrete
#recorder Element -file $dataDir/MVLEM_panel_stress_steel.out -time -ele 1 RCPanel 1 panel_stress_steel

# Unaxial Steel Recorders for all panels
#recorder Element -file $dataDir/MVLEM_strain_stress_steel1_1.out -time -ele 1 RCPanel 1 strain_stress_steelX
#recorder Element -file $dataDir/MVLEM_strain_stress_steel2_1.out -time -ele 1 RCPanel 1 strain_stress_steelY
#recorder Element -file $dataDir/MVLEM_strain_stress_steel1_2.out -time -ele 1 RCPanel 2 strain_stress_steelX
#recorder Element -file $dataDir/MVLEM_strain_stress_steel2_2.out -time -ele 1 RCPanel 2 strain_stress_steelY
#recorder Element -file $dataDir/MVLEM_strain_stress_steel1_3.out -time -ele 1 RCPanel 3 strain_stress_steelX
#recorder Element -file $dataDir/MVLEM_strain_stress_steel2_3.out -time -ele 1 RCPanel 3 strain_stress_steelY
#recorder Element -file $dataDir/MVLEM_strain_stress_steel1_4.out -time -ele 1 RCPanel 4 strain_stress_steelX
#recorder Element -file $dataDir/MVLEM_strain_stress_steel2_4.out -time -ele 1 RCPanel 4 strain_stress_steelY
#recorder Element -file $dataDir/MVLEM_strain_stress_steel1_5.out -time -ele 1 RCPanel 5 strain_stress_steelX
#recorder Element -file $dataDir/MVLEM_strain_stress_steel2_5.out -time -ele 1 RCPanel 5 strain_stress_steelY

#recorder Element -file $dataDir/MVLEM_strain_stress_concrete1.out -time -ele 1 RCPanel 1 strain_stress_concrete1
#recorder Element -file $dataDir/MVLEM_strain_stress_concrete2.out -time -ele 1 RCPanel 1 strain_stress_concrete2
#recorder Element -file $dataDir/MVLEM_strain_stress_interlock1.out -time -ele 1 RCPanel 1 strain_stress_interlock1
#recorder Element -file $dataDir/MVLEM_strain_stress_interlock2.out -time -ele 1 RCPanel 1 strain_stress_interlock2

# Cracking angles for all panels
#recorder Element -file $dataDir/MVLEM_cracking_angles1.out -time -ele 1 2 3 4 5 RCPanel 1 cracking_angles
#recorder Element -file $dataDir/MVLEM_cracking_angles2.out -time -ele 1 2 3 4 5 RCPanel 2 cracking_angles
#recorder Element -file $dataDir/MVLEM_cracking_angles3.out -time -ele 1 2 3 4 5 RCPanel 3 cracking_angles
#recorder Element -file $dataDir/MVLEM_cracking_angles4.out -time -ele 1 2 3 4 5 RCPanel 4 cracking_angles
#recorder Element -file $dataDir/MVLEM_cracking_angles5.out -time -ele 1 2 3 4 5 RCPanel 5 cracking_angles

# ---------------------
# Define Axial Load
# ---------------------

set N [expr 149.0]; # kips

# -------------------------------------------------------
# Set parameters for displacement controlled analysis
# -------------------------------------------------------

# vector of displacement-cycle peaks in terms of wall drift ratio (TOTAL displacements)
set iDmax "0.001 0.0025 0.005 0.0075 0.01 0.015 0.02 0.03"; 
set Ncycles 1;	# specify the number of cycles at each peak
set CycleType Full;		# type of static analysis: Full / Push / Half
set TotalDrift [expr 4*$Ncycles*(0.001 + 0.0025 + 0.005 + 0.0075 + 0.01 + 0.015 + 0.02 + 0.03)]

#set iDmax "0.030";
#set Ncycles 1;	# specify the number of cycles at each peak
#set CycleType Push;		# type of static analysis: Full / Push / Half
#set TotalDrift [expr 0.030]

set TotalSteps [expr $TotalDrift*$H/$Dincr]
	
set Tol 1.0e-5;
set LunitTXT "inch";

puts "Model generated sucessfully"