# Guanzhou modified to test 27node

# Simple Shear loading a 27-node brick element
# [gjie,jeremic]@ucdavis.edu
# 31Oct2003
# for solid models, 27 node brick



# ################################
# create the modelbuilder
# #################################

model BasicBuilder -ndm 3 -ndf 3

set p1   0.05555556
set np1 -0.05555556
set p2   0.11111111
set np2 -0.11111111
set p4   0.44444444
set np4 -0.44444444
set p3   0.22222222
set np3 -0.22222222

set g   9.81
set lf1 1.0

# ################################
# build the model
# #################################

node  1 1.0 0.707106781 0.707106781
node  4 1.0 0.0 0.0
node  5 1.0 1.414213562 0.0
node  8 1.0 0.707106781 -0.707106781
node  9 1.0 1.060660172 0.35355339
node 12 1.0 0.35355339 -0.35355339
node 16 1.0 0.35355339 0.35355339
node 20 1.0 1.060660172 -0.35355339
node 24 1.0 0.707106781 0.0

node 13 0.5 0.707106781 0.707106781
node 15 0.5 0.0 0.0
node 17 0.5 1.414213562 0.0
node 19 0.5 0.707106781 -0.707106781
node 21 0.5 1.060660172 0.35355339
node 23 0.5 0.35355339 -0.35355339
node 25 0.5 0.35355339 0.35355339
node 26 0.5 1.060660172 -0.35355339
node 27 0.5 0.707106781 0.0

node  2 0.0 0.707106781 0.707106781
node  3 0.0 0.0 0.0
node  6 0.0 1.414213562 0.0
node  7 0.0 0.707106781 -0.707106781
node 10 0.0 1.060660172 0.35355339
node 11 0.0 0.35355339 -0.35355339
node 14 0.0 0.35355339 0.35355339
node 18 0.0 1.060660172 -0.35355339
node 22 0.0 0.707106781 0.0

fix  3 1 1 1
fix  4 0 1 1
fix 15 0 1 1
fix  5 0 0 1
fix 17 0 0 1
fix  6 1 0 1
fix  2 1 0 0
fix 14 1 0 0
fix 10 1 0 0
fix 18 1 0 0
fix  7 1 0 0
fix 11 1 0 0
fix 22 1 0 0

# elastic material
nDMaterial ElasticIsotropic3D 1 70000 0.0 1.8

#(28 args)______tag________________________27 nodes______________________________________matID_bforce1_bforce2_bforce3_Rho
element Brick27N 1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 1    0.0     0.0    -9.81   1.8

#===========================================================
# Lateral loading

pattern Plain 2 Linear {
   load  1 0 0 $np1 
   load  2 0 0 $np1 
   load  5 0 $p1 0
   load  6 0 $p1 0
   load  7 0 0 $p1
   load  8 0 0 $p1
   load  9 0 $p2 $np2
   load 10 0 $p2 $np2
   load 11 0 $np2 $p2
   load 12 0 $np2 $p2
   load 13 0 0 $np3
   load 14 0 $np2 $np2
   load 16 0 $np2 $np2
   load 17 0 $p3 0
   load 18 0 $p2 $p2
   load 19 0 0 $p3
   load 20 0 $p2 $p2
   load 21 0 $p4 $np4
   load 23 0 $np4 $p4
   load 25 0 $np4 $np4
   load 26 0 $p4 $p4
}

# ----------------------------
# Start of recorder generation
# ----------------------------

#recorder display ShakingBeam 0 0 300 300 -wipe
#prp -100 100 120.5
#vup 0 1 1 
#display 1 0 1 


# ################################
# create the analysis
# #################################

system UmfPack
constraints Penalty 1e12 1e12
#test NormDispIncr 1.0e-8 30 1
test NormUnbalance 1.0e-10 30 1
integrator LoadControl $lf1 1 $lf1 $lf1
algorithm Newton
numberer RCM
analysis Static

for {set i 1} {$i <=1} {incr i} {
 puts $i
 analyze 1
}


print node 1 5 8

wipe
