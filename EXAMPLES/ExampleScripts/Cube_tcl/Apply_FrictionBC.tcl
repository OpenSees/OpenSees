#================================================================
# Apply_FrictionBC.tcl
# -- Apply frition B.C.'s at top and bottom platen
# Zhaohui Yang UC Davis
# April 6, 2002
#================================================================

#Remove Isotropic consolidation B.C.'s prior to applying friction
for {set i 0} {$i < $SPC_counter} {incr i} {
  remove SPconstraint $i
}
puts "Finished removing old $SPC_counter SPconstraints..."

#Open file nodes.dat again to apply frition constrains
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
  
  if { $Z == 0 } {
    fix  $tag 1 1 1
  }
  if { $Z == $Zmax } {
    fix  $tag 1 1 0
  }  
}
puts " Finished applying friction to the top and bottom platens..."
close $nodes

