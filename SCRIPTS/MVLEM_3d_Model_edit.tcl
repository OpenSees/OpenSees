# ----------------------------------------------------------------
# Example 1: Simulation of wall flexural response using MVLEM
# Specimen: RW2 (Thomsen and Wallace, 1995)
# Created by: Kristijan Kolozvari (kkolozvari@fullerton.edu)
# Date: 7/2015
# ----------------------------------------------------------------

# --------------------------------------------------------
# Start of model generation (Units: kip, in, sec, ksi)
# --------------------------------------------------------
wipe;
# Set Up Directories
#set modelName "RW2"; 			# Model Name
#set dataDir MVLEM_$modelName;	# Name of output folder
#file mkdir $dataDir;

# Create ModelBuilder for 2D element (with two-dimensions and 2 DOF/node)
model BasicBuilder -ndm 3 -ndf 6

# --------------------------------------------
# Set geometry, nodes, boundary conditions
# --------------------------------------------

# Wall Geometry 
set H 144.0; 	# Wall height
set t 4.0;		# Wall thickness

# Create nodes 
# node nodeId xCrd yCrd

node 1	[expr 0*48.0/1]	[expr 0*$H/8]	0.0;  
node 2	[expr 1*48.0/1]	[expr 0*$H/8]	0.0;

node 11	[expr 0*48.0/1]	[expr 1*$H/8]	0.0;  
node 12	[expr 1*48.0/1]	[expr 1*$H/8]	0.0;

node 21	[expr 0*48.0/1]	[expr 2*$H/8]	0.0;  
node 22	[expr 1*48.0/1]	[expr 2*$H/8]	0.0;

node 31	[expr 0*48.0/1]	[expr 3*$H/8]	0.0;  
node 32	[expr 1*48.0/1]	[expr 3*$H/8]	0.0;

node 41	[expr 0*48.0/1]	[expr 4*$H/8]	0.0;  
node 42	[expr 1*48.0/1]	[expr 4*$H/8]	0.0;

node 51	[expr 0*48.0/1]	[expr 5*$H/8]	0.0;  
node 52	[expr 1*48.0/1]	[expr 5*$H/8]	0.0;

node 61	[expr 0*48.0/1]	[expr 6*$H/8]	0.0;  
node 62	[expr 1*48.0/1]	[expr 6*$H/8]	0.0;

node 71	[expr 0*48.0/1]	[expr 7*$H/8]	0.0;  
node 72	[expr 1*48.0/1]	[expr 7*$H/8]	0.0;

node 81	[expr 0*48.0/1]	[expr 8*$H/8]	0.0;  
node 82	[expr 1*48.0/1]	[expr 8*$H/8]	0.0;

# Boundary conditions
fix 1 1 1 1 1 1 1;	
fix 2 1 1 1 1 1 1;	

# Set Control Node and DOF
set IDctrlNode 82;
set IDctrlDOF 1;

# ------------------------------------------------------------------------
# Define uniaxial materials 
# ------------------------------------------------------------------------

# STEEL ...........................................................
# uniaxialMaterial SteelMPF $mattag  $fyp  $fyn  $E0  $bp  $bn   $R0  $a1  $a2

# steel Y boundary
set fyYbp 57.3; 
set bybp 0.0185; 
set fyYbn 63.0; 
set bybn 0.02; 		

# steel Y web
set fyYwp 48.8; 	
set bywp 0.035; 
set fyYwn 65.0; 	
set bywn 0.02; 		

# steel misc
set Es 29000.0; 	
set R0 20.0; 		
set a1 0.925; 		
set a2 0.0015;		

# Build steel materials
uniaxialMaterial    SteelMPF  1 $fyYbp $fyYbn $Es $bybp $bybn $R0 $a1 $a2;
uniaxialMaterial    SteelMPF  2 $fyYwp $fyYwn $Es $bywp $bywn $R0 $a1 $a2;

# Set MVLEM Reinforcing Ratios
set rouYb 0.029333; 
set rouYw 0.003333; 

# CONCRETE ........................................................
# uniaxialMaterial ConcreteCM $mattag  $fpcc  $epcc  $Ec  $rc  $xcrn   $ft  $et  $rt   $xcrp <-GapClose $gap>

# unconfined
set fpc 6.2; 	
set ec0 -0.0021;	
set ft 0.295;		
set et 0.00008;		
set Ec 4500; 		
set xcrnu 1.039;	
set xcrp 10000;		
set ru 7;			
set rt 1.2;			

# confined
set fpcc 6.9036; 	
set ec0c -0.0033;	
set Ecc 5091.3; 	
set xcrnc 1.0125;	
set rc 7.3049;		

# Build concrete materials
# confined concrete
uniaxialMaterial ConcreteCM 3 -$fpcc  $ec0c  $Ecc $rc  $xcrnc  $ft  $et  $rt  $xcrp -GapClose 1; 
# unconfined concrete
uniaxialMaterial ConcreteCM 4 -$fpc   $ec0   $Ec  $ru  $xcrnu  $ft  $et  $rt  $xcrp -GapClose 1; 


# CONCRETE ........................................................
#uniaxialMaterial Concrete02 $matTag $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets
# unconfined
# set fpc 6.2; 		# peak compressive stress
# set ec0 -0.0021;	# strain at peak compressive stress
# set ft 0.295;		# peak tensile stress
# set et 0.00008;		# strain at peak tensile stress

# confined
# set fpcc 6.9036; 	# peak compressive stress
# set ec0c -0.0033;	# strain at peak compressive stress

# Bould concrete materials
# uniaxialMaterial Concrete02 4 -$fpc   $ec0  [expr -0.1*$fpc] -0.012 0.3 $ft 100; # unconfined concrete
# uniaxialMaterial Concrete02 3 -$fpcc   $ec0c  [expr -0.1*$fpcc] -0.080 0.3 $ft 100; # confined concrete

# SHEAR ........................................................
# uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
# NOTE: large shear stiffness assigned since only flexural response
set Ashweb [expr 4.0*30.0];					
set G [expr 0.1*0.4*57.0* pow($fpc*1000.0,0.5)];	
set GAs [expr $G * $Ashweb]; 	

# Build shear material
uniaxialMaterial Elastic 5 $GAs;

# ------------------------------
#  Define MVLEM elements
# ------------------------------

#element FourNodeMVLEM3D eleTag iNode jNode   kNode lNode   m      c     nu Tfactor -thick fiberThick              -width fiberWidth                      -rho Rho                                               -matConcrete matTagsConcrete -matSteel matTagsSteel    -matShear matTagShear\n";
element FourNodeMVLEM3D	11	0.0	1	2	12	11	8	0.4	0.2	0.6	-thick $t $t $t $t $t $t $t $t -width 7.5 1.5 7.5 7.5 7.5 7.5 1.5 7.5 -rho $rouYb 0.0 $rouYw $rouYw $rouYw $rouYw 0.0 $rouYb -matConcrete 3 4 4 4 4 4 4 3 -matSteel 1 2 2 2 2 2 2 1 -matShear 5

element FourNodeMVLEM3D	21	0.0	11	12	22	21	8	0.4	0.2	0.6	-thick $t $t $t $t $t $t $t $t -width 7.5 1.5 7.5 7.5 7.5 7.5 1.5 7.5 -rho $rouYb 0.0 $rouYw $rouYw $rouYw $rouYw 0.0 $rouYb -matConcrete 3 4 4 4 4 4 4 3 -matSteel 1 2 2 2 2 2 2 1 -matShear 5

element FourNodeMVLEM3D	31	0.0	21	22	32	31	8	0.4	0.2	0.6	-thick $t $t $t $t $t $t $t $t -width 7.5 1.5 7.5 7.5 7.5 7.5 1.5 7.5 -rho $rouYb 0.0 $rouYw $rouYw $rouYw $rouYw 0.0 $rouYb -matConcrete 3 4 4 4 4 4 4 3 -matSteel 1 2 2 2 2 2 2 1 -matShear 5

element FourNodeMVLEM3D	41	0.0	31	32	42	41	8	0.4	0.2	0.6	-thick $t $t $t $t $t $t $t $t -width 7.5 1.5 7.5 7.5 7.5 7.5 1.5 7.5 -rho $rouYb 0.0 $rouYw $rouYw $rouYw $rouYw 0.0 $rouYb -matConcrete 3 4 4 4 4 4 4 3 -matSteel 1 2 2 2 2 2 2 1 -matShear 5

element FourNodeMVLEM3D	51	0.0	41	42	52	51	8	0.4	0.2	0.6	-thick $t $t $t $t $t $t $t $t -width 7.5 1.5 7.5 7.5 7.5 7.5 1.5 7.5 -rho $rouYb 0.0 $rouYw $rouYw $rouYw $rouYw 0.0 $rouYb -matConcrete 3 4 4 4 4 4 4 3 -matSteel 1 2 2 2 2 2 2 1 -matShear 5

element FourNodeMVLEM3D	61	0.0	51	52	62	61	8	0.4	0.2	0.6	-thick $t $t $t $t $t $t $t $t -width 7.5 1.5 7.5 7.5 7.5 7.5 1.5 7.5 -rho $rouYb 0.0 $rouYw $rouYw $rouYw $rouYw 0.0 $rouYb -matConcrete 3 4 4 4 4 4 4 3 -matSteel 1 2 2 2 2 2 2 1 -matShear 5

element FourNodeMVLEM3D	71	0.0	61	62	72	71	8	0.4	0.2	0.6	-thick $t $t $t $t $t $t $t $t -width 7.5 1.5 7.5 7.5 7.5 7.5 1.5 7.5 -rho $rouYb 0.0 $rouYw $rouYw $rouYw $rouYw 0.0 $rouYb -matConcrete 3 4 4 4 4 4 4 3 -matSteel 1 2 2 2 2 2 2 1 -matShear 5

element FourNodeMVLEM3D	81	0.0	71	72	82	81	8	0.4	0.2	0.6	-thick $t $t $t $t $t $t $t $t -width 7.5 1.5 7.5 7.5 7.5 7.5 1.5 7.5 -rho $rouYb 0.0 $rouYw $rouYw $rouYw $rouYw 0.0 $rouYb -matConcrete 3 4 4 4 4 4 4 3 -matSteel 1 2 2 2 2 2 2 1 -matShear 5

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
#recorder Node -file $dataDir/MVLEM_DOFs.out -time -node 1 2 3 4 -dof 1 2 3 disp

# Element recorders
#recorder Element -file MVLEM_Fgl.out -time -ele 11 globalForce
#recorder Element -file MVLEM_Curvature.out -time -ele 11 Curvature

# Fiber responses
#recorder Element -file MVLEM_fiber_strain.out -time -ele 11 fiber_strain 
#recorder Element -file MVLEM_fiber_stress_concrete.out -time -ele 11 fiber_stress_concrete
#recorder Element -file MVLEM_fiber_stress_steel.out -time -ele 11 fiber_stress_steel

# Shear spring response
#recorder Element -file MVLEM_shear_force_def.out -time -ele 11 shear_force_deformation

#recorder Node -file MVLEM_base.out -time -node 1 2 -dof 1 2 3 4 5 6 reaction

# ---------------------
# Define Axial Load
# ---------------------

set N [expr 85.0];
print node $IDctrlNode
# -------------------------------------------------------
# Set parameters for displacement controlled analysis
# -------------------------------------------------------

# vector of displacement-cycle peaks in terms of wall drift ratio (flexural displacements)
#set iDmax "0.000330792 0.001104233 0.002925758 0.004558709 0.006625238 0.010816268 0.014985823 0.019655056";
set iDmax "0.02";
set Dincr 0.02;		
set CycleType Push;		
set Ncycles 1;			
set Tol 1.0e-5;
set LunitTXT "inch";