model basic -ndm 3 -ndf 6

# create the material
section ElasticMembranePlateSection 1 3.0e3 0.25 1.175 1.27

# set some parameters for node and element generation
set Plate ShellMITC4

set eleArgs "1"

#these should both be even
set nx 8
set ny 2

#loaded nodes
set mid [expr ( ($nx+1)*($ny+1)+1 ) / 2 ]
set side1 [expr ($nx + 2)/2 ]
set side2 [expr ($nx+1)*($ny+1) - $side1 + 1 ]

# generate the nodes and elements
block2D $nx $ny 1 1 $Plate $eleArgs {
    1 -20 0 0
    2 -20 0 40
    3 20 0 40
    4 20 0 0
    5 -10 10 20
    7 10 10 20
    9 0 10 20
}

# Set some parameters
set outFile ARL360.g3
source ReadSMDFile.tcl
ReadSMDFile ARL360.at2 $outFile dt
timeSeries Path 1 -filePath $outFile -dt $dt -factor 1.0
pattern UniformExcitation  2   1  -accel 1


# define the boundary conditions
# rotation free about x-axis (remember right-hand-rule)
fixZ 0.0 1 1 1 0 1 1
fixZ 40.0 1 1 1 0 1 1

system BandGeneral
constraints Plain
test NormUnbalance 1.0e-8  10 1
algorithm Newton
numberer RCM
integrator Newmark  0.5  0.25 
analysis Transient

set tFinal [expr 20000 * $dt]
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
puts "eigen values at start of transient: [eigen -Umfpack 1]"
puts "eigen values at start of transient: [eigen -Umfpack 2]"
puts "eigen values at start of transient: [eigen -Umfpack 2]"

# Print state of node 3
#print node 3
#print ele 1

puts [getTime]
wipe

