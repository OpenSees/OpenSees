#================================================================
# InputElements.tcl
# --Create all elements from file elements.dat
# Zhaohui Yang and Boris Jeremic (@ucdavis.edu)
# March 31, 2002
#================================================================

#File to store elements which belong to pile
set els_in_pile [open "els_in_pile.dat" "w"]
puts $els_in_pile 46; # I know there are 46 elememts in the pile

#Open file elements.dat containing element info
set elements [open "elements.dat" "r"]
gets $elements num_of_element
#puts " Number of element is $num_of_element ."
					      
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
  set n9   [lindex $line  9];
  set n10  [lindex $line 10];
  set n11  [lindex $line 11];
  set n12  [lindex $line 12];
  set n13  [lindex $line 13];
  set n14  [lindex $line 14];
  set n15  [lindex $line 15];
  set n16  [lindex $line 16];
  set n17  [lindex $line 17];
  set n18  [lindex $line 18];
  set n19  [lindex $line 19];
  set n20  [lindex $line 20];
  set gar  [lindex $line 21];
  set Mat_no [lindex $line 22];

  #Create the element
  #_________________tag__________________________20__...__nodes_________________________________________________matID_bforce1,2&3___rho
  if {($Mat_no) == 2 } {
     element Brick20N $tag $n5 $n6 $n7 $n8 $n1 $n2 $n3 $n4 $n13 $n14 $n15 $n16 $n9 $n10 $n11 $n12 $n17 $n18 $n19 $n20 3 0.0 0.0 $g $rho_pile
  } elseif { ($Mat_no) == 3 } {
     element Brick20N $tag $n5 $n6 $n7 $n8 $n1 $n2 $n3 $n4 $n13 $n14 $n15 $n16 $n9 $n10 $n11 $n12 $n17 $n18 $n19 $n20 3 0.0 0.0 $g 0.0
  } else {
     element Brick20N $tag $n5 $n6 $n7 $n8 $n1 $n2 $n3 $n4 $n13 $n14 $n15 $n16 $n9 $n10 $n11 $n12 $n17 $n18 $n19 $n20 1 0.0 0.0 $g $rho_sand
  }
  
  #If this is a pile element, save it in els_in_pile.dat, for replacing after self-weight loading
  if { ($Mat_no) == 2 } {
    puts $els_in_pile "$tag $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 $n9 $n10 $n11 $n12 $n13 $n14 $n15 $n16 $n17 $n18 $n19 $n20 $Mat_no"
  }
}
puts "Finished creating $num_of_element elements..."
close $els_in_pile
close $elements


