model basic -ndm 3 -ndf 6
set PI 3.14159
set g 386.4

set fc 6.0
set massD [expr 0.15/(12.*12.*12*$g)]
set E [expr 57000*sqrt($fc)]

# create the material
section ElasticMembranePlateSection 1 $E 0.25 12.0 $massD

# set some parameters for node and element generation
set Plate ShellMITC4

set eleArgs "1"

#these should both be even
set nx 8
set ny 2


#loaded nodes
set top [expr ( ($nx+1)*($ny+1))]

# generate the nodes and elements
block2D $nx $ny 1 1 $Plate $eleArgs {
    1 -120 0 0
    2 -120 0 360
    3  120 0 360
    4  120 0 0
}

# Set some parameters
set outFile elCentro.g3
source ReadSMDFile.tcl
ReadSMDFile elCentro.AT2 $outFile dt
timeSeries Path 1 -filePath $outFile -dt $dt -factor $g
pattern UniformExcitation  2   1  -accel 1

# define the boundary conditions
# rotation free about x-axis (remember right-hand-rule)
fixZ 0.0 1 1 1 1 1 1
#fixZ 360. 1 1 1 0 1 1

recorder Node -file node.out -time -node $top -dof 1 disp

system BandGeneral
constraints Plain
test NormUnbalance 1.0e-8  10 0
algorithm Newton
numberer RCM
integrator Newmark  0.5  0.25 
analysis Transient

set tFinal [expr 1559 * $dt]
set tCurrent [getTime]
set ok 0

# Perform the transient analysis
while {$ok == 0 && $tCurrent < $tFinal} {
    
    set ok [analyze 1 .01]
    
    # if the analysis fails try initial tangent iteration
    if {$ok != 0} {
	puts "regular newton failed .. lets try an initail stiffness for this step"
	test NormDispIncr 1.0e-12  100 0
	algorithm ModifiedNewton -initial
	set ok [analyze 1 .01]
	if {$ok == 0} {puts "that worked .. back to regular newton"}
	test NormDispIncr 1.0e-12  10 
	algorithm Newton
    }
    
    set tCurrent [getTime]
}

# Print a message to indicate if analysis successful or not
if {$ok == 0} {
   puts "Transient analysis completed SUCCESSFULLY";
} else {
   puts "Transient analysis completed FAILED";    
}

# Perform an eigenvalue analysis
set lambda [eigen 1]
puts " Period T1: [expr 2*$PI/sqrt([lindex $lambda 0])] sec (frequency: [lindex $lambda 0])"    

# Print state of node 3
#print node 3
#print ele 1

puts [getTime]
wipe

