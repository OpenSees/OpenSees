#================================================================
# InputNodes.tcl
# --Create the all nodes and SP_Constraints from file nodes.dat
# Zhaohui Yang UCDavis
# April 6, 2002
#================================================================

#Open file nodes.dat containing node info
set nodes [open "nodes.dat" "r"]
gets $nodes num_of_node
					      
#Open file nodes.ops containing node info for visualizing
set onodes [open "nodes.ops" "w"]
puts $onodes $num_of_node

#Initialize sp constraint counters
set SPC_counter	0;

set Zmax 0; # Var to store the max. Z coor.
for {set i 0} {$i < $num_of_node} {incr i} {
  gets $nodes line
  set tag  [lindex $line  0];
  set Tx   [lindex $line  1];
  set Ty   [lindex $line  2];
  set Tz   [lindex $line  3];
  set Rx   [lindex $line  4];
  set Ry   [lindex $line  5];
  set Rz   [lindex $line  6];
  set Xs   [lindex $line  7];  set X   [expr $Xs*$m ];
  set Ys   [lindex $line  8];  set Y   [expr $Ys*$m ];
  set Zs   [lindex $line  9];  set Z   [expr $Zs*$m ];
  
  #Create the node
  node $tag $X $Y $Z
  
  #Output node info to nodes.ops
  puts $onodes "$tag $Tx $Ty $Tz $Rx $Ry $Rz $X $Y $Z"

  #Create the SP_Constraints for the above node
  set temp [expr $Tx + $Ty + $Tz]
  #puts "Sum of Txyz = $temp."
  if { $temp > 0 } {
    fix  $tag $Tx $Ty $Tz

    #Count how many sp constraint created
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
  
  #Get the max. Z
  if { abs($Z) >  abs($Zmax) } {
    set Zmax $Z
  }
}
puts " Finished creating $num_of_node nodes..."
puts " Finished creating $SPC_counter constraints..."
close $nodes
#Check Zmax
puts " Zmax = $Zmax"

