

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
ops.model('BasicBuilder','-ndm',3,'-ndf',6)

# --------------------------------------------
# Set geometry, nodes, boundary conditions
# --------------------------------------------

# Wall Geometry 
ops.set('H','144.0;','#','Wall','height')
ops.set('t','4.0;','#','Wall','thickness')

# Create nodes 
# node nodeId xCrd yCrd

ops.node(1,'[expr','0*48.0/1]','[expr','0*$H/8]','0.0;')
ops.node(2,'[expr','1*48.0/1]','[expr','0*$H/8]','0.0;')

ops.node(11,'[expr','0*48.0/1]','[expr','1*$H/8]','0.0;')
ops.node(12,'[expr','1*48.0/1]','[expr','1*$H/8]','0.0;')

ops.node(21,'[expr','0*48.0/1]','[expr','2*$H/8]','0.0;')
ops.node(22,'[expr','1*48.0/1]','[expr','2*$H/8]','0.0;')

ops.node(31,'[expr','0*48.0/1]','[expr','3*$H/8]','0.0;')
ops.node(32,'[expr','1*48.0/1]','[expr','3*$H/8]','0.0;')

ops.node(41,'[expr','0*48.0/1]','[expr','4*$H/8]','0.0;')
ops.node(42,'[expr','1*48.0/1]','[expr','4*$H/8]','0.0;')

ops.node(51,'[expr','0*48.0/1]','[expr','5*$H/8]','0.0;')
ops.node(52,'[expr','1*48.0/1]','[expr','5*$H/8]','0.0;')

ops.node(61,'[expr','0*48.0/1]','[expr','6*$H/8]','0.0;')
ops.node(62,'[expr','1*48.0/1]','[expr','6*$H/8]','0.0;')

ops.node(71,'[expr','0*48.0/1]','[expr','7*$H/8]','0.0;')
ops.node(72,'[expr','1*48.0/1]','[expr','7*$H/8]','0.0;')

ops.node(81,'[expr','0*48.0/1]','[expr','8*$H/8]','0.0;')
ops.node(82,'[expr','1*48.0/1]','[expr','8*$H/8]','0.0;')

# Boundary conditions
ops.fix(1,1,1,1,1,1,'1;','#','Fixed','condition','at','node',1)
ops.fix(2,1,1,1,1,1,'1;','#','Fixed','condition','at','node',2)

# Set Control Node and DOF
ops.set('IDctrlNode','82;')
ops.set('IDctrlDOF','1;')

# ------------------------------------------------------------------------
# Define uniaxial materials 
# ------------------------------------------------------------------------

# STEEL ...........................................................
# uniaxialMaterial SteelMPF $mattag  $fyp  $fyn  $E0  $bp  $bn   $R0  $a1  $a2

# steel Y boundary
ops.set('fyYbp','57.3;','#','fy','-','tension')
ops.set('bybp','0.0185;','#','strain','hardening','-','tension')
ops.set('fyYbn','63.0;','#','fy','-','compression')
ops.set('bybn','0.02;','#','strain','hardening','-','compression')

# steel Y web
ops.set('fyYwp','48.8;','#','fy','-','tension')
ops.set('bywp','0.035;','#','strain','hardening','-','tension')
ops.set('fyYwn','65.0;','#','fy','-','compression')
ops.set('bywn','0.02;','#','strain','hardening','-','compression')

# steel misc
ops.set('Es','29000.0;','#','Young's','modulus')
ops.set('R0','20.0;','#','initial','value','of','curvature','parameter')
ops.set('a1','0.925;','#','curvature','degradation','parameter')
ops.set('a2','0.0015;','#','curvature','degradation','parameter')

# Build steel materials
ops.uniaxialMaterial('SteelMPF',1,'$fyYbp','$fyYbn','$Es','$bybp','$bybn','$R0','$a1','$a2;','#','steel','Y','boundary')
ops.uniaxialMaterial('SteelMPF',2,'$fyYwp','$fyYwn','$Es','$bywp','$bywn','$R0','$a1','$a2;','#','steel','Y','web')

# Set MVLEM Reinforcing Ratios
ops.set('rouYb','0.029333;','#','Y','boundary')
ops.set('rouYw','0.003333;','#','Y','web')

# CONCRETE ........................................................
# uniaxialMaterial ConcreteCM $mattag  $fpcc  $epcc  $Ec  $rc  $xcrn   $ft  $et  $rt   $xcrp <-GapClose $gap>

# unconfined
ops.set('fpc','6.2;','#','peak','compressive','stress')
ops.set('ec0','-0.0021;','#','strain','at','peak','compressive','stress')
ops.set('ft','0.295;','#','peak','tensile','stress')
ops.set('et','0.00008;','#','strain','at','peak','tensile','stress')
ops.set('Ec','4500;','#','Young's','modulus')
ops.set('xcrnu','1.039;','#','cracking','strain','-','compression')
ops.set('xcrp','10000;','#','cracking','strain','-','tension')
ops.set('ru','7;','#','shape','parameter','-','compression')
ops.set('rt','1.2;','#','shape','parameter','-','tension')

# confined
ops.set('fpcc','6.9036;','#','peak','compressive','stress')
ops.set('ec0c','-0.0033;','#','strain','at','peak','compressive','stress')
ops.set('Ecc','5091.3;','#','Young's','modulus')
ops.set('xcrnc','1.0125;','#','cracking','strain','-','compression')
ops.set('rc','7.3049;','#','shape','parameter','-','compression')

# Build concrete materials
# confined concrete
ops.uniaxialMaterial('ConcreteCM',3,'-$fpcc','$ec0c','$Ecc','$rc','$xcrnc','$ft','$et','$rt','$xcrp','-GapClose','1;')
# unconfined concrete
ops.uniaxialMaterial('ConcreteCM',4,'-$fpc','$ec0','$Ec','$ru','$xcrnu','$ft','$et','$rt','$xcrp','-GapClose','1;')


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

# Combine concrete and steel
ops.nDMaterial('FSAM',101,0.0,2,1,3,0.0050,'$rouYb',0.2,'0.01;')
ops.nDMaterial('FSAM',102,0.0,2,1,3,0.0025,0.0,0.2,'0.01;')
ops.nDMaterial('FSAM',103,0.0,2,2,4,0.0025,'$rouYw',0.2,'0.01;')
# ------------------------------
#  Define SFI elements
# ------------------------------

ops.element('FourNodeSFI_MVLEM3D',11,1,2,12,11,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-mat',101,102,103,103,103,103,102,101)
ops.element('FourNodeSFI_MVLEM3D',21,11,12,22,21,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-mat',101,102,103,103,103,103,102,101)
ops.element('FourNodeSFI_MVLEM3D',31,21,22,32,31,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-mat',101,102,103,103,103,103,102,101)
ops.element('FourNodeSFI_MVLEM3D',41,31,32,42,41,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-mat',101,102,103,103,103,103,102,101)
ops.element('FourNodeSFI_MVLEM3D',51,41,42,52,51,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-mat',101,102,103,103,103,103,102,101)
ops.element('FourNodeSFI_MVLEM3D',61,51,52,62,61,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-mat',101,102,103,103,103,103,102,101)
ops.element('FourNodeSFI_MVLEM3D',71,61,62,72,71,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-mat',101,102,103,103,103,103,102,101)
ops.element('FourNodeSFI_MVLEM3D',81,71,72,82,81,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-mat',101,102,103,103,103,103,102,101)

# ------------------------------
# End of model generation
# ------------------------------

# Initialize
initialize

# ------------------------------
# Recorder generation
# ------------------------------

# Nodal recorders
ops.recorder('Node','-file','SFI_MVLEM_Dtop.out','-time','-node','$IDctrlNode','-dof',1,'disp')
#recorder Node -file $dataDir/SFI_MVLEM_DOFs.out -time -node 1 2 3 4 -dof 1 2 3 disp

# Element recorders
#recorder Element -file SFI_MVLEM_Fgl.out -time -ele 11 globalForce
#recorder Element -file SFI_MVLEM_Curvature.out -time -ele 11 Curvature

# Fiber responses
#recorder Element -file SFI_MVLEM_fiber_strain_L.out -time -ele 11 RCPanel 1 panel_strain
#recorder Element -file SFI_MVLEM_fiber_strain_R.out -time -ele 11 RCPanel 8 panel_strain
# recorder Element -file SFI_MVLEM_fiber_strain.out -time -ele 11 fiber_strain 
# recorder Element -file SFI_MVLEM_fiber_stress_concrete.out -time -ele 11 fiber_stress_concrete
# recorder Element -file SFI_MVLEM_fiber_stress_steel.out -time -ele 11 fiber_stress_steel

# Shear spring response
# recorder Element -file SFI_MVLEM_shear_force_def.out -time -ele 11 shear_force_deformation

#recorder Node -file SFI_MVLEM_base.out -time -node 1 2 -dof 1 2 3 4 5 6 reaction

# ---------------------
# Define Axial Load
# ---------------------

ops.set('N','[expr','85.0];','#','kips')
ops.printModel('node','$IDctrlNode')
# -------------------------------------------------------
# Set parameters for displacement controlled analysis
# -------------------------------------------------------

# vector of displacement-cycle peaks in terms of wall drift ratio (flexural displacements)
#set iDmax "0.000330792 0.001104233 0.002925758 0.004558709 0.006625238 0.010816268 0.014985823 0.019655056";
ops.set('iDmax','"6.405093125000000e-04',0.001615589090278,0.003556868736111,0.005507881097222,0.007798747881944,0.012445314194444,0.008062608909722,0.012470400916667,0.016982636347222,'0.021999289805556";')
ops.set('Dincr','0.031;','#','displacement','increment','for','displacement','controlled','analysis.')
ops.set('CycleType','Full;','#','type','of','static','analysis:','Full','/','Push','/','Half')
ops.set('Ncycles','2;','#','specify','the','number','of','cycles','at','each','peak')
ops.set('Tol','1.0e-7;')
ops.set('LunitTXT','"inch";')
