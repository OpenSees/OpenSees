# Lateral loading of a 20-node brick element
# zhyang@ucdavis.edu
# jeremic@ucdavis.edu
# 24Aug2001
# for solid models, 20 node brick



# ################################
# create the modelbuilder
# #################################

model BasicBuilder -ndm 3 -ndf 3

set p   1.667
set np -6.667
set g   9.81
set lf1 1.0

# ################################
# build the model
# #################################

node  1 1.0 1.0 0.0
node  2 0.0 1.0 0.0
node  3 0.0 0.0 0.0
node  4 1.0 0.0 0.0
node  5 1.0 1.0 2.0
node  6 0.0 1.0 2.0
node  7 0.0 0.0 2.0
node  8 1.0 0.0 2.0
node  9 0.5 1.0 0.0
node 10 0.0 0.5 0.0
node 11 0.5 0.0 0.0
node 12 1.0 0.5 0.0
node 13 0.5 1.0 2.0
node 14 0.0 0.5 2.0
node 15 0.5 0.0 2.0
node 16 1.0 0.5 2.0
node 17 1.0 1.0 1.0
node 18 0.0 1.0 1.0
node 19 0.0 0.0 1.0
node 20 1.0 0.0 1.0

fix  1 1 1 1
fix  2 1 1 1  
fix  3 1 1 1
fix  4 1 1 1
fix  9 1 1 1
fix 10 1 1 1
fix 11 1 1 1
fix 12 1 1 1
    
# elastic material
nDMaterial ElasticIsotropic3D 1 70000 0.3 1.8

#(28 args)______tag________________________20 nodes___________________matID_bforce1_bforce2_bforce3_Rho
element Brick20N 1  5 6 7 8 1 2 3 4 13 14 15 16 9 10 11 12 17 18 19 20  1     0.0     0.0    -9.81  1.8

#===========================================================
# Lateral loading

pattern Plain 2 Linear {
   load  5 0 $p  0 
   load  6 0 $p  0 
   load  7 0 $p  0 
   load  8 0 $p  0 
   load 13 0 $np 0 
   load 14 0 $np 0 
   load 15 0 $np 0 
   load 16 0 $np 0 
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


print node 5 6 7 8

wipe
