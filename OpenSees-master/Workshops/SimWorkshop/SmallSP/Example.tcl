model basic -ndm 3 -ndf 3
nDMaterial ElasticIsotropic   1   100   0.25  1.27

set eleArgs "1" 
set element stdBrick

set nz 10
set nx 4
set ny 4

set nn [expr ($nz+1)*($nx+1)*($ny+1) ]

# mesh generation
block3D $nx $ny $nz   1 1  $element  $eleArgs {
    1   -1     -1      0
    2    1     -1      0
    3    1      1      0
    4   -1      1      0 
    5   -1     -1     10
    6    1     -1     10
    7    1      1     10
    8   -1      1     10
}


set load 0.10
pattern Plain 1 Linear {
   load $nn  $load  $load  0.0
}

fixZ 0.0   1 1 1 

integrator LoadControl  1.0  1 
test NormUnbalance     1.0e-10    20     0
algorithm Newton
numberer RCM
constraints Plain 
system Mumps
analysis Static 

analyze 5

recorder Node -file Node.out -time -node $nn -dof 1 disp

# Remove the static analysis & reset the time to 0.0
wipeAnalysis
setTime 0.0

remove loadPattern 1


rayleigh 0.01 0.0 0.0 0.0

# Create the transient analysis
test EnergyIncr     1.0e-10    20   0
algorithm Newton
numberer RCM
constraints Plain 
integrator Newmark 0.5 0.25
system Mumps
analysis Transient


analyze 1000 1.0

puts [nodeDisp $nn]



