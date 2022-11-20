#================================================================
# Plug_In_PileElements.tcl
# --Replace soil elements in pile position with real pile elements
# Zhaohui Yang and Boris Jeremic (@ucdavis.edu)
# March 31, 2002
#================================================================

#File to store elements which belong to pile
set els_in_pile [open "els_in_pile.dat" "r"]
gets $els_in_pile num_of_pile_element;
					      
for {set i 0} {$i < $num_of_pile_element} {incr i} {
  gets $els_in_pile line
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
  set Mat_no [lindex $line 21];

  #Remove dummy soil elements first
  remove element $tag 
  #puts "  removed $tag dummy element..."

  #Plug in pile elements
  #_________________tag__________________________20__...__nodes_________________________________________________matID_bforce1,2&3___rho
  if {$Mat_no == 2 } {
     element Brick20N $tag $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 $n9 $n10 $n11 $n12 $n13 $n14 $n15 $n16 $n17 $n18 $n19 $n20 3 0.0 0.0 $g $rho_pile
     #puts "  added $tag pile element..."
  }   

}
puts "Finished Plugging in pile elements..."
close $els_in_pile
