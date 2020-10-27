wipe
model basic -ndm 3 -ndf 3
# MATERIALS
set fcc		-30

#nDMaterial Damage2p 1	$fcc -fct $fct -H $H -theta $theta -rho_bar $rho_bar -E $E -Gt $Gt -Gc $Gc -ni $ni -tangent 0
nDMaterial Damage2p 1	$fcc

#nDMaterial ElasticIsotropic 1 $Ec  $ni 

set nx 1
set ny 1
set nz 1
set e1 1
set n1 1
block3D $nx $ny $nz $e1 $n1 stdBrick 1 {
  1    0    0    0
  2 1000    0    0
  3 1000 1000    0
  4    0 1000    0
  5    0    0 1000
  6 1000    0 1000
  7 1000 1000 1000
  8    0 1000 1000
}

#printGID check.msh

fix 1	1 1 1
fix 3	1 0 1
fix 5	1 0 0
fix 7	1 0 0

#Recorders
recorder Node -file damage2p.txt -time -node 8 -dof 1 disp

#Loads
set V 0.25e6
pattern Plain 1 Linear {
	load 2  $V  0.0  0.0
	load 4  $V  0.0  0.0
	load 6  $V  0.0  0.0
	load 8  $V  0.0  0.0
}

set startTime [clock clicks -milliseconds]

system SparseGeneral -piv
test EnergyIncr 1e-16 5 0
numberer RCM
constraints Plain
algorithm Newton


integrator DisplacementControl 8 1 0.0 ;
analysis Static
#initialize
#analyze 18
#analyze 40
analyze 1

integrator DisplacementControl 8 1 -0.2
analyze 40
integrator DisplacementControl 8 1 0.2
analyze 30
integrator DisplacementControl 8 1 -0.2
analyze 50


set finishTime [clock clicks -milliseconds];
set timeSeconds  [expr ($finishTime-$startTime)/1000];
set timeMinutes  [expr ($timeSeconds/60)];
set timeHours    [expr ($timeSeconds/3600)];
set timeMinutes  [expr ($timeMinutes - $timeHours*60)];
set timeSeconds  [expr ($timeSeconds - $timeMinutes*60 - $timeHours*3600)];
puts "----------------------------------";
puts "TOTAL TIME TAKEN $timeHours:$timeMinutes:$timeSeconds";
puts "\n";



