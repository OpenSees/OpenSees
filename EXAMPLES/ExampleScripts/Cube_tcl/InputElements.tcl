#================================================================
# InputElements.tcl
# --Create all elements from file elements.dat
# Zhaohui Yang
# April 6, 2002
#================================================================

#Open file elements.dat containing element info
set elements [open "elements.dat" "r"]
gets $elements num_of_element
puts " Number of element is $num_of_element ."
					      
#Open file elements.ops containing element info for visualizing
set oelements [open "elements.ops" "w"]
puts $oelements "$num_of_element $nonper"

for {set i 0} {$i < $num_of_element} {incr i} {
  gets $elements line
  set tag  [lindex $line  0];
  set n1   [lindex $line  1];
  set n2   [lindex $line  2];
  set n3   [lindex $line  3];
  set n4   [lindex $line  4];
  set n5   [lindex $line  5];
  set n6   [lindex $line  6];
  set n7   [lindex $line  7];
  set n8   [lindex $line  8];
  set gar  [lindex $line 9];
  set Mat_no [lindex $line 10];

  #puts "Adding element $tag, $n1, $n4, $n5, $n8, $Mat_no..."
  #Create the element
  #_________________tag__________________________20__...__nodes_________________________________________________matID_bforce1,2&3___rho
  element Brick8N $tag $n5 $n6 $n7 $n8 $n1 $n2 $n3 $n4 1 0.0 0.0 $g $rho_sand
  #element Brick8N $tag $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 1 0.0 0.0 $g $rho_sand

  #Output elements for visualizing
  puts $oelements "$tag $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 $gar $Mat_no"

}
puts " Finished creating $num_of_element elements..."
close $elements

