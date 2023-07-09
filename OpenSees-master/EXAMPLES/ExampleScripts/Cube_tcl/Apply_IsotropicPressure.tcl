#================================================================
# Apply_IsotropicPressure.tcl
# --Create isotropic load on all surface nodes from  surfnodes.dat
# Zhaohui Yang UCDavis
# April 6, 2002
#================================================================

#stage 1--isotropic load
#===========================================================
#Open file containing surface nodes for applying confining pressure
set surfnds [open "surfnodes.dat" "r"]

#Assign isotropic confining pressure which will be applied on the sample
set P [expr 100.0*$kPa]

#Read in number of node on surface & number of unit force on top surface
gets $surfnds line
set num_surfnds [lindex $line 0];
set Nforce      [lindex $line 1];

set  p [ expr $P / $Nforce ]

puts " Number of nodes on surface is $num_surfnds ."
puts " Unit force is $p ."

# Create nodal loads at nodes for isotropic consolidation
pattern Plain 1 Linear {
   # Read the surface loads and apply them
   for {set i 0} {$i < $num_surfnds} {incr i} {
     gets $surfnds line
     set tag  [lindex $line  0];
     set n1   [lindex $line  1];
     set n2   [lindex $line  2];
     set n3   [lindex $line  3];
     #puts " loading Node $tag with Fx = $n1 Fy = $n2 Fz = $n3"
     #Create the nodal load
     #______tag_________Fx___________Fy_____________Fz
     load  $tag  [expr $n1*$p] [expr $n2*$p] [expr $n3*$p]
   }

}

puts " Finished applying loads on $num_surfnds surface nodes..."
close $surfnds

