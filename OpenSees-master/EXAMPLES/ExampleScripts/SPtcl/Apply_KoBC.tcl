#================================================================
# Apply_KoBC.tcl
# --Create Ko consolidation B.C.'s
# Zhaohui Yang and Boris Jeremic (@ucdavis.edu)
# March 31, 2002
#================================================================

#Open file nodes.dat containing node info
set nodes [open "nodes.dat" "r"]	      
gets $nodes num_of_node
set SPC_counter	0;

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
  
  # Fix bottom nodes and leave other nodes only move vertically
  if { $Z == $Zbot } {
    set Tx 1;
    set Ty 1;
    set Tz 1;
  } else {
    set Tx 1;
    set Ty 1;
    set Tz 0;
  }

  #Create the SP_Constraints for the above node
  fix  $tag [expr $Tx] [expr $Ty] [expr $Tz]  
  if { $Tx } {
     set SPC_counter [expr $SPC_counter+1]
  }
  if { $Ty } {
     set SPC_counter [expr $SPC_counter+1]
  }
  if { $Tz } {
     set SPC_counter [expr $SPC_counter+1]
  }

}

puts "Finished creating $SPC_counter SPconstraints for Ko consolidation B.C.'s..."
close $nodes

