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
ReadRecord $record.AT2  $record.dat dT nPts	    

timeSeries Path 2 -filePath $record.dat -dt $dT 
pattern UniformExcitation  2 1 -accel 2

rayleigh 0.0 0.0 0.0 0.0

wipeAnalysis

system BandGeneral
constraints Plain
test NormDispIncr 1.0e-8  10 
algorithm Newton
numberer RCM
integrator Newmark  0.5  0.25 
analysis Transient

# create a recorder
recorder Node -time -file disp.out -node 3 4 -dof 1 2 3 disp
recorder Node -time -file a2.out -timeSeries 1 -node 3 4 -dof 1 2 3 accel
recorder Node -time -file a1.out -node 3 4 -dof 1 2 3 accel

set tFinal [expr $nPts * $dT]
set tCurrent [getTime]
set ok 0

while {$ok == 0 && $tCurrent < $tFinal} {
    
    set ok [analyze 1 $dT]
    
    # if the analysis fails try initial tangent iteration
    if {$ok != 0} {
	puts "regular newton failed .. lets try an initail stiffness for this step"
	test NormDispIncr 1.0e-8  1000 1
	algorithm ModifiedNewton -initial
	set ok [analyze 1 $dT]
	if {$ok == 0} {puts "that worked .. back to regular newton"}
	test NormDispIncr 1.0e-12  10 
	algorithm Newton
    }
    
    set tCurrent [getTime]
}

# Print a message to indicate if analysis succesfull or not
if {$ok == 0} {
   puts "Transient analysis completed SUCCESSFULLY";
} else {
   puts "Transient analysis completed FAILED";    
}
