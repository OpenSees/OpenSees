# Do operations of Example3.1 by sourcing in the tcl file
source RCFrameGravity.tcl
puts "Gravity Analysis Completed"

# Set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0

# Define Pattern for Lateral Reference loads
set H 10.0;  
pattern Plain 2 1 {
     load 3 $H 0.0 0.0 
     load 4 $H 0.0 0.0 
}

set dU 0.1;  
integrator DisplacementControl  3   1   $dU  1 $dU $dU

recorder Node -file node32.out -time -node 3 4 -dof 1 2 3 disp
recorder EnvelopeElement -file ele12.out -time -ele 1 2 forces

# Set some parameters
set maxU 15.0;	        # Max displacement
set currentDisp 0.0;
set ok 0

while {$ok == 0 && $currentDisp < $maxU} {
	set ok [analyze 1]

	# if the analysis fails try initial tangent iteration
	if {$ok != 0} {
	    puts "regular newton failed .. lets try an initial stiffness for this step"
	    test NormDispIncr 1.0e-12  1000 
	    algorithm ModifiedNewton -initial
	    set ok [analyze 1]
	    test NormDispIncr 1.0e-12  10 
	    algorithm Newton
	}

	set currentDisp [nodeDisp 3 1]
}


if {$ok == 0} {
  puts "Pushover analysis completed SUCCESSFULLY";
} else {
  puts "Pushover analysis FAILED";    
}

# Print the state at node 3
print node 3
