#[jeremic,zhyang]@ucdavis.edu
# tcl example file using Template3Dep material
# -- Manzari-Dafalias material model -- loose of critical sands
# cf. Geotechnique v.47 No.2 255-272, 1997
# Inelastic static analysis
#
#two load stages:
# stage 1: isotropic compression
# stage 2: triaxial shearing
# Feb. 24, 2003


# ################################
# create the modelbuilder
# #################################

model BasicBuilder -ndm 3 -ndf 3

set g  -9.81
set p   25.0
set np -25.0
set lf1 0.05
set NN1 10
# ################################
# build the model
# #################################

node 1 1.0 1.0 0.0
node 2 0.0 1.0 0.0
node 3 0.0 0.0 0.0
node 4 1.0 0.0 0.0
node 5 1.0 1.0 1.0
node 6 0.0 1.0 1.0
node 7 0.0 0.0 1.0
node 8 1.0 0.0 1.0

#triaxial test boundary
fix 1 0 1 1
fix 2 1 1 1
fix 3 1 0 1
fix 4 0 0 1
 

#Elastic-plastic model
#================================================================
# Manzari-Dafalias model
#     var   Manzari-Dafalias Yield Surface
set DPys "-MD"
#     var   Manzari-Dafalias Potential Surfac
set DPps "-MD"
#     var Scalar hardening law for m
set ES1  "-Leq  0.0"
#     var bounding_surface hardening rule for back-stress alpha_ij
set ET1  "-NLEijMD -Mc 1.14 -Me 1.14 -kc_b 3.975 -ke_b 4.20 -ke_b 2.0 -ke_d 0.07 -ho 1200 -Cm 0.0 -Ao 2.64 -Fmax 100 -Cf 100"
# initial stress
set stressp "-0.10 0 0  0 -0.10 0  0 0 -0.10"

#____________E______Eo_____v__rho________________po
set EPS "70000.0 70000.0 0.3 1.8 -NOD 3 -NOS 2 0.05 0.0 -stressp $stressp -e 0.8 -ec 0.75 -lambda 0.025 -p_ref 100.0 -a 0.5"

#
nDMaterial Template3Dep 1 -YS $DPys -PS $DPps -EPS $EPS -ELS1 $ES1 -ELT1 $ET1

#____________tag_____8 nodes______matID__bforce1_bforce2_bforce3_rho
element Brick8N  1  5 6 7 8 1 2 3 4   1      0.0     0.0    $g   1.8
    
#stage 1
#===========================================================
#isotropic load
pattern Plain 2 "Linear" {
load 1 $np 0   0   -pattern 2
load 3 0   $p  0   -pattern 2
load 4 $np $p  0   -pattern 2
load 5 $np $np $np -pattern 2
load 6 $p  $np $np -pattern 2
load 7 $p  $p  $np -pattern 2
load 8 $np $p  $np -pattern 2
}

##Set up recorder
#recorder Node node_iso.out disp -time -node 5 -dof 3
#recorder plot node_iso.out "disp load"  10 10 300 300 -columns 2 1

#stress recorder
recorder Element all -time -file SigmaIso.out stress
#recorder plot SigmaIso.out "Load factor vs. Normal stress"  10 200 300 300 -columns 1 3

system UmfPack
constraints Penalty 1e12 1e12
test NormDispIncr 1.0e-8 30 0
integrator LoadControl $lf1 1 $lf1 $lf1
algorithm Newton
numberer RCM
analysis Static

analyze $NN1


#stage 2
#===========================================================
#Axial loading
wipeAnalysis
#loadConst

equalDOF 5 6 3 3
equalDOF 5 7 3 3
equalDOF 5 8 3 3
       	  

#set previous load constant
loadConst time 0

# Before shifting to disp. control, apply constant P shearing
set p   -50.0
set mp   50.0
set mq -100.0

pattern Plain 3 "Linear" {
load 1 $mp 0.0 0.0 
load 3 0.0 $p  0.0 
load 4 $mp $p  0.0 
load 5 $mp $mp $mq
load 6 $p  $mp $mq
load 7 $p  $p  $mq 
load 8 $mp $p  $mq
}

# Stop the old recorders by destroying them
remove recorders

#Set up recorder for axial loading stage
#recorder Node node_z.out disp -time -node 5 -dof 3

recorder Node node_z.out disp -time -node 5 -dof 3
#recorder plot node_z.out "Disp vs Pseudo-Load" 10 100 300 300 -columns 2 1

recorder Element 1 -file element.out pqall
#recorder plot element.out "p vs. q" 10 500 300 300 -columns 1 2

set NN2 300
set ndz -0.0002
system UmfPack
constraints Penalty 1e12 1e12
test NormDispIncr 1.0e-08 30  0
integrator DisplacementControl 5 3 $ndz 1 $ndz $ndz
algorithm Newton
numberer RCM
analysis Static

analyze $NN2

wipe
