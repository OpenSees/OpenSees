# Do operations of Example3.1 by sourcing in the tcl file
source RCFrameGravity.tcl
puts "Gravity load analysis completed"

# Set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0

# Define nodal mass in terms of axial load on columns
set g 386.4
set m [expr ($W/2.0)/$g];       # expr command to evaluate an expression

#    tag   MX   MY   RZ
mass  3    $m   $m    1.0e-16
mass  4    $m   $m    1.0e-16

# Define dynamic loads
# --------------------

# Set some parameters
set record IELC180

# Source in TCL proc to read PEER SMD record
source ReadRecord.tcl
ReadRecord $record.DT2  $record.dat dT nPts	    

timeSeries Path 2 -filePath $record.dat -dt $dT -factor $cm
pattern MultiSupport  2   {
    groundMotion 5 Plain -disp 2
    imposedMotion 1 1 5
    imposedMotion 2 1 5
}

# set the rayleigh damping factors for nodes & elements
rayleigh 0.0 0.0 0.0 0.0

# remove sp constraints from nodes where motion to be input
remove sp 1 1
remove sp 2 1

wipeAnalysis

system BandGeneral
constraints Transformation
test NormDispIncr 1.0e-8  10 
algorithm Newton
numberer RCM
integrator Newmark  0.5  0.25 
analysis Transient

recorder Node -time -file multi.out -node 1 3 -dof 1 disp

set tFinal [expr $nPts * $dT]
set tCurrent [getTime]
set ok 0


while {$ok == 0 && $tCurrent < $tFinal} {
    set ok [analyze 1 $dT]
    
    # if the analysis fails try initial tangent iteration
    if {$ok != 0} {
	puts "regular newton failed .. lets try an initial stiffness for this step"
	test NormDispIncr 1.0e-8  1000 1
	algorithm ModifiedNewton -initial
	set ok [analyze 1 $dT]
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

