# ----------------------------------------------------------------
# Example 1: Simulation of wall flexural response using MVLEM
# Specimen: RW2 (Thomsen and Wallace, 1995)
# Created by: Kristijan Kolozvari (kkolozvari@fullerton.edu)
# Date: 7/2015
# ----------------------------------------------------------------

# --------------------------------------------------------
# Start of model generation (Units: kip, in, sec, ksi)
# --------------------------------------------------------
wipe
# Set Up Directories
# modelName "RW2" 			# Model Name
# dataDir MVLEM_$modelName	# Name of output folder
#file mkdir $dataDir

# Create ModelBuilder for 2D element (with two-dimensions and 2 DOF/node)
model('BasicBuilder','-ndm',3,'-ndf',6)

# --------------------------------------------
# Set geometry, nodes, boundary conditions
# --------------------------------------------

# Wall Geometry 
H = 144.0 #Wallheight
t = 4.0 #Wallthickness

# Create nodes 
# node nodeId xCrd yCrd

node(1,0,0,0.0)
node(2,48,0,0.0)

node(11,0,18,0.0)
node(12,48,18,0.0)

node(21,0,36,0.0)
node(22,48,36,0.0)

node(31,0,54,0.0)
node(32,48,54,0.0)

node(41,0,72,0.0)
node(42,48,72,0.0)

node(51,0,90,0.0)
node(52,48,90,0.0)

node(61,0,108,0.0)
node(62,48,108,0.0)

node(71,0,126,0.0)
node(72,48,126,0.0)

node(81,0,144,0.0)
node(82,48,144,0.0)

# Boundary conditions
fix(1,1,1,1,1,1,1) #Fixedconditionatnode1)
fix(2,1,1,1,1,1,1) #Fixedconditionatnode2)

# Set Control Node and DOF
IDctrlNode = 82)
IDctrlDOF = 1)

# ------------------------------------------------------------------------
# Define uniaxial materials 
# ------------------------------------------------------------------------

# STEEL ...........................................................
# uniaxialMaterial SteelMPF $mattag  $fyp  $fyn  $E0  $bp  $bn   $R0  $a1  $a2

# steel Y boundary
fyYbp = 57.3 #fy-tension
bybp = 0.0185 #strainhardening-tension
fyYbn = 63.0 #fy-compression
bybn = 0.02 #strainhardening-compression

# steel Y web
fyYwp = 48.8 #fy-tension
bywp = 0.035 #strainhardening-tension
fyYwn = 65.0 #fy-compression
bywn = 0.02 #strainhardening-compression

# steel misc
Es = 29000.0 #Young'smodulus
R0 = 20.0 #initialvalueofcurvatureparameter
a1 = 0.925 #curvaturedegradationparameter
a2 = 0.0015 #curvaturedegradationparameter

# Build steel materials
uniaxialMaterial('SteelMPF',1,'$fyYbp','$fyYbn','$Es','$bybp','$bybn','$R0','$a1','$a2') #steelYboundary)
uniaxialMaterial('SteelMPF',2,'$fyYwp','$fyYwn','$Es','$bywp','$bywn','$R0','$a1','$a2') #steelYweb)

# Set MVLEM Reinforcing Ratios
rouYb = 0.029333 #Yboundary
rouYw = 0.003333 #Yweb

# CONCRETE ........................................................
# uniaxialMaterial ConcreteCM $mattag  $fpcc  $epcc  $Ec  $rc  $xcrn   $ft  $et  $rt   $xcrp <-GapClose $gap>

# unconfined
fpc = 6.2 #peakcompressivestress
ec0 = -0.0021 #strainatpeakcompressivestress
ft = 0.295 #peaktensilestress
et = 0.00008 #strainatpeaktensilestress
Ec = 4500 #Young'smodulus
xcrnu = 1.039 #crackingstrain-compression
xcrp = 10000 #crackingstrain-tension
ru = 7 #shapeparameter-compression
rt = 1.2 #shapeparameter-tension

# confined
fpcc = 6.9036 #peakcompressivestress
ec0c = -0.0033 #strainatpeakcompressivestress
Ecc = 5091.3 #Young'smodulus
xcrnc = 1.0125 #crackingstrain-compression
rc = 7.3049 #shapeparameter-compression

# Build concrete materials
# confined concrete
uniaxialMaterial('ConcreteCM',3,'-$fpcc','$ec0c','$Ecc','$rc','$xcrnc','$ft','$et','$rt','$xcrp','-GapClose',1)
# unconfined concrete
uniaxialMaterial('ConcreteCM',4,'-$fpc','$ec0','$Ec','$ru','$xcrnu','$ft','$et','$rt','$xcrp','-GapClose',1)


# CONCRETE ........................................................
#uniaxialMaterial Concrete02 $matTag $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets
# unconfined
#  fpc 6.2 		# peak compressive stress
#  ec0 -0.0021	# strain at peak compressive stress
#  ft 0.295		# peak tensile stress
#  et 0.00008		# strain at peak tensile stress

# confined
#  fpcc 6.9036 	# peak compressive stress
#  ec0c -0.0033	# strain at peak compressive stress

# Bould concrete materials
# uniaxialMaterial Concrete02 4 -$fpc   $ec0  [expr -0.1*$fpc] -0.012 0.3 $ft 100 # unconfined concrete
# uniaxialMaterial Concrete02 3 -$fpcc   $ec0c  [expr -0.1*$fpcc] -0.080 0.3 $ft 100 # confined concrete

# SHEAR ........................................................
# uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
# NOTE: large shear stiffness assigned since only flexural response
Ashweb = 120 #Grossareaofthewallcrosssection
G = 179 #ShearModulus0.1G
GAs = 120*179 #ShearStiffness

# Build shear material
uniaxialMaterial('Elastic',5,'$GAs')

# ------------------------------
#  Define MVLEM elements
# ------------------------------

#element FourNodeMVLEM3D eleTag iNode jNode   kNode lNode   m      c     nu Tfactor -thick fiberThick              -width fiberWidth                      -rho Rho                                               -matConcrete matTagsConcrete -matSteel matTagsSteel    -matShear matTagShear\n"
element('FourNodeMVLEM3D',11,0.0,1,2,12,11,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-rho','$rouYb',0.0,'$rouYw','$rouYw','$rouYw','$rouYw',0.0,'$rouYb','-matConcrete',3,4,4,4,4,4,4,3,'-matSteel',1,2,2,2,2,2,2,1,'-matShear',5)

element('FourNodeMVLEM3D',21,0.0,11,12,22,21,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-rho','$rouYb',0.0,'$rouYw','$rouYw','$rouYw','$rouYw',0.0,'$rouYb','-matConcrete',3,4,4,4,4,4,4,3,'-matSteel',1,2,2,2,2,2,2,1,'-matShear',5)

element('FourNodeMVLEM3D',31,0.0,21,22,32,31,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-rho','$rouYb',0.0,'$rouYw','$rouYw','$rouYw','$rouYw',0.0,'$rouYb','-matConcrete',3,4,4,4,4,4,4,3,'-matSteel',1,2,2,2,2,2,2,1,'-matShear',5)

element('FourNodeMVLEM3D',41,0.0,31,32,42,41,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-rho','$rouYb',0.0,'$rouYw','$rouYw','$rouYw','$rouYw',0.0,'$rouYb','-matConcrete',3,4,4,4,4,4,4,3,'-matSteel',1,2,2,2,2,2,2,1,'-matShear',5)

element('FourNodeMVLEM3D',51,0.0,41,42,52,51,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-rho','$rouYb',0.0,'$rouYw','$rouYw','$rouYw','$rouYw',0.0,'$rouYb','-matConcrete',3,4,4,4,4,4,4,3,'-matSteel',1,2,2,2,2,2,2,1,'-matShear',5)

element('FourNodeMVLEM3D',61,0.0,51,52,62,61,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-rho','$rouYb',0.0,'$rouYw','$rouYw','$rouYw','$rouYw',0.0,'$rouYb','-matConcrete',3,4,4,4,4,4,4,3,'-matSteel',1,2,2,2,2,2,2,1,'-matShear',5)

element('FourNodeMVLEM3D',71,0.0,61,62,72,71,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-rho','$rouYb',0.0,'$rouYw','$rouYw','$rouYw','$rouYw',0.0,'$rouYb','-matConcrete',3,4,4,4,4,4,4,3,'-matSteel',1,2,2,2,2,2,2,1,'-matShear',5)

element('FourNodeMVLEM3D',81,0.0,71,72,82,81,8,0.4,0.2,0.6,'-thick','$t','$t','$t','$t','$t','$t','$t','$t','-width',7.5,1.5,7.5,7.5,7.5,7.5,1.5,7.5,'-rho','$rouYb',0.0,'$rouYw','$rouYw','$rouYw','$rouYw',0.0,'$rouYb','-matConcrete',3,4,4,4,4,4,4,3,'-matSteel',1,2,2,2,2,2,2,1,'-matShear',5)

# ------------------------------
# End of model generation
# ------------------------------

# Initialize
initialize

# ------------------------------
# Recorder generation
# ------------------------------

# Nodal recorders
recorder('Node','-file','MVLEM_Dtop.out','-time','-node','$IDctrlNode','-dof',1,'disp')
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

N = 85.0 #kips
printModel('node','$IDctrlNode')
# -------------------------------------------------------
# Set parameters for displacement controlled analysis
# -------------------------------------------------------

# vector of displacement-cycle peaks in terms of wall drift ratio (flexural displacements)
# iDmax "0.000330792 0.001104233 0.002925758 0.004558709 0.006625238 0.010816268 0.014985823 0.019655056"
iDmax = '"0.02"')
Dincr = 0.02 #displacementincrementfordisplacementcontrolledanalysis.
CycleType = 'Push' #typeofstaticanalysis:Full/Push/Half
Ncycles = 1 #specifythenumberofcyclesateachpeak
Tol = 1.0e-5)
LunitTXT = '"inch"')
