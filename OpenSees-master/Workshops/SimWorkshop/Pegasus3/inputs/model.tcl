#UNITS
  set in 1.0
  set kip 1.0
  set g 386.4
#NODES
 node 1 0.0 0.0 -mass 0.0 0.0 0.0
 node 2 360.0 0.0 -mass 0.0 0.0 0.0
 node 3 720.0 0.0 -mass 0.0 0.0 0.0
 node 4 1080.0 0.0 -mass 0.0 0.0 0.0
 node 5 0.0 180.0 -mass 0.42054865424430643 0.42054865424430643 0.0
 node 6 360.0 180.0 -mass 0.42054865424430643 0.42054865424430643 0.0
 node 7 720.0 180.0 -mass 0.42054865424430643 0.42054865424430643 0.0
 node 8 1080.0 180.0 -mass 0.42054865424430643 0.42054865424430643 0.0
 node 9 0.0 360.0 -mass 0.42054865424430643 0.42054865424430643 0.0
 node 10 360.0 360.0 -mass 0.42054865424430643 0.42054865424430643 0.0
 node 11 720.0 360.0 -mass 0.42054865424430643 0.42054865424430643 0.0
 node 12 1080.0 360.0 -mass 0.42054865424430643 0.42054865424430643 0.0
 node 13 0.0 540.0 -mass 0.2911490683229814 0.2911490683229814 0.0
 node 14 360.0 540.0 -mass 0.2911490683229814 0.2911490683229814 0.0
 node 15 720.0 540.0 -mass 0.2911490683229814 0.2911490683229814 0.0
 node 16 1080.0 540.0 -mass 0.2911490683229814 0.2911490683229814 0.0
source SteelWSections.tcl
#SECTIONS
SteelWSection 1 W14x370 1
SteelWSection 4 W33x141 1
SteelWSection 2 W14x370 1
SteelWSection 5 W33x130 1
SteelWSection 3 W14x211 1
SteelWSection 6 W27x102 1
#Transformations
geomTransf Linear 1
geomTransf PDelta 2
geomTransf Corotational 3
#Column Elements
  element dispBeamColumn 1 1 5 5 1 1
  element dispBeamColumn 2 5 9 5 2 1
  element dispBeamColumn 3 9 13 5 3 1
  element dispBeamColumn 4 2 6 5 1 1
  element dispBeamColumn 5 6 10 5 2 1
  element dispBeamColumn 6 10 14 5 3 1
  element dispBeamColumn 7 3 7 5 1 1
  element dispBeamColumn 8 7 11 5 2 1
  element dispBeamColumn 9 11 15 5 3 1
  element dispBeamColumn 10 4 8 5 1 1
  element dispBeamColumn 11 8 12 5 2 1
  element dispBeamColumn 12 12 16 5 3 1
#Beam Elements
  element dispBeamColumn 13 5 6 5 4 1
  element dispBeamColumn 14 6 7 5 4 1
  element dispBeamColumn 15 7 8 5 4 1
  element dispBeamColumn 16 9 10 5 5 1
  element dispBeamColumn 17 10 11 5 5 1
  element dispBeamColumn 18 11 12 5 5 1
  element dispBeamColumn 19 13 14 5 6 1
  element dispBeamColumn 20 14 15 5 6 1
  element dispBeamColumn 21 15 16 5 6 1
#Constraints
  fix 1 1 1 1
  fix 2 1 1 1
  fix 3 1 1 1
  fix 4 1 1 1
#Gravity Loads
  timeSeries Linear 1
  pattern Plain 1 1 {
    load 1 0.0 -0.0 0.0
    load 2 0.0 -0.0 0.0
    load 3 0.0 -0.0 0.0
    load 4 0.0 -0.0 0.0
    load 5 0.0 -162.5 0.0
    load 6 0.0 -162.5 0.0
    load 7 0.0 -162.5 0.0
    load 8 0.0 -162.5 0.0
    load 9 0.0 -162.5 0.0
    load 10 0.0 -162.5 0.0
    load 11 0.0 -162.5 0.0
    load 12 0.0 -162.5 0.0
    load 13 0.0 -112.5 0.0
    load 14 0.0 -112.5 0.0
    load 15 0.0 -112.5 0.0
    load 16 0.0 -112.5 0.0
}
# GRAVITY ANALYSIS
  numberer RCM; 
  test NormDispIncr 1.0e-12 6 0
  algorithm Newton; 
  constraints Plain
  system UmfPack; 
  integrator LoadControl 0.1;
  analysis Static
  set ok [analyze 10]
# load const & return nodes
  loadConst -time 0.0
set node "16 3 3"
