#================================================================
# InputNodes.tcl
# --Create the all nodes and SP_Constraints from file nodes.dat
# Zhaohui Yang UC Davis
# March 31, 2002
#================================================================

#Open file nodes.dat containing node info
set nodes [open "$Dir/nodes.dat" "r"]	      
gets $nodes num_of_node

#Open file nodes.ops containing node info for visualizing
set onodes [open "$Dir/nodes.ops" "w"]
puts $onodes $num_of_node
					      
set Zbot 0; # Var to store the deepest Z coor.
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
  
  #Create the node
  #Mesh is in meter, so specify with [ *$m ]. If in inch, please specify with [ *$in ]
  node $tag [expr $X*$m ] [expr $Y*$m ] [expr $Z*$m ]

  #Output node info to nodes.ops
  puts $onodes "$tag $Tx $Ty $Tz $Rx $Ry $Rz $X $Y $Z"
  
  #Find the deepest Z
  if { abs($Z) >  abs($Zbot) } {
    set Zbot $Z
  }
}
puts "Finished creating $num_of_node nodes..."
close $nodes
close $onodes
#Check Zbot
#puts "Zbot = $Zbot"

