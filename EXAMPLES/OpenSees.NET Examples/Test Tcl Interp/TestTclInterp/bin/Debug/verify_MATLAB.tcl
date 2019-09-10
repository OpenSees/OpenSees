
model basic -ndm 2
set PI 3.1415
set nd1 1
node $nd1 0. 0. 0
fix $nd1 0 1 1

set nd0 2
node $nd0 0. -1 0
fix $nd0 1 1 1

set uy [expr 10./1000];
set m  [expr 50.96e6/200]
set K1 [expr 0.92e7]
set K2 0.92e6;
set Alpha [expr $K2/$K1];
set Fy 92000.;

mass $nd1 $m $m 0

uniaxialMaterial Elastic 1 $K1
set eta 100
set beta 0.5
set gamma 0.5
# element elastomericBearingBoucWen 1 $nd0 $nd1 $K1 $Fy $Alpha 0 0 $eta $beta $gamma -P 1 -Mz 1
uniaxialMaterial Steel01 2 $Fy $K1 $Alpha
element twoNodeLink 1 $nd0 $nd1 -mat 1 2 1 -dir 1 2 3 
timeSeries Path 1 -dt 0.02 -filePath accel.txt -factor 9.81
pattern UniformExcitation  1 1 -accel 1

recorder Node -file disp.txt -time -node $nd1 -dof 1 disp
recorder Node -file react.txt -time -node $nd0 -dof 1 reaction
system BandGeneral
constraints Transformation
numberer RCM 
test NormDispIncr 1e-8 100 0
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient
analyze 2000 0.02


