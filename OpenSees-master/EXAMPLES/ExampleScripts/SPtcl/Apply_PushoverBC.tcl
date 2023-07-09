#================================================================
# Apply_PushoverBC.tcl
# --Apply Pushover B.C.'s
# Zhaohui Yang and Boris Jeremic (@ucdavis.edu)
# March 31, 2002
#================================================================


#Remove Ko B.C.'s prior to applying Pushover B.C.'s

for {set i 0} {$i < $SPC_counter} {incr i} {
  remove SPconstraint $i
}
puts "Finished removing old $SPC_counter SPconstraints..."


#Open file nodes.dat containing node info
set nodes [open "nodes.dat" "r"]	      
gets $nodes num_of_node
					      
for {set i 0} {$i < $num_of_node} {incr i} {
  gets $nodes line
  set tag  [lindex $line  0];
  set Tx   [lindex $line  1];
  set Ty   [lindex $line  2];
  set Tz   [lindex $line  3];
  set Rx   [lindex $line  4];
  set Ry   [lindex $line  5];
  set Rz   [lindex $line  6];
  set Xs   [lindex $line  7];  set X   [expr $Xs  ];
  set Ys   [lindex $line  8];  set Y   [expr $Ys  ];
  set Zs   [lindex $line  9];  set Z   [expr $Zs  ];
  

  #Create the SP_Constraints for the above node
  fix  $tag [expr $Tx] [expr $Ty] [expr $Tz]  
}

puts "Finished creating Pushover B.C.'s..."
close $nodes

