#================================================================
# InputElements.tcl
# --Create all elements from file elements.dat
# Zhaohui Yang  UC Davis
# March 31, 2002
#================================================================


#Material number for plastic bowl elements
set PBMat 2

#File to store elements which belong to the plastic bowl
set els_in_pb [open "$Dir/PBElements.dat" "w"]
puts $els_in_pb 89; # I know there are 89 elememts in the plastic bowl

#Open file elements.dat containing element info
set elements [open "$Dir/elements.dat" "r"]
gets $elements num_of_element
#puts " Number of element is $num_of_element ."

#Open file elements.ops containing element info for visualizing
set oelements [open "$Dir/elements.ops" "w"]
puts $oelements "$num_of_element $nonper"
					      
for {set i 0} {$i < $num_of_element} {incr i} {
  gets $elements line
  if { $nonper == 20} {
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
  } elseif { $nonper == 8 } {
    set tag  [lindex $line  0];
    set n1   [lindex $line  1];
    set n2   [lindex $line  2];
    set n3   [lindex $line  3];
    set n4   [lindex $line  4];
    set n5   [lindex $line  5];
    set n6   [lindex $line  6];
    set n7   [lindex $line  7];
    set n8   [lindex $line  8];
    set gar  [lindex $line  9];
    set Mat_no [lindex $line 10];
  } else {
    puts "Currently can only handle 8 or 20 node brick element."
    exit;
  }
      
  #Create the element
  #_________________tag__________________________20__...__nodes_________________________________________________matID_bforce1,2&3___rho
  if { $nonper == 20} {  
    if {($Mat_no) == 1 } {
       element Brick20N $tag $n5 $n6 $n7 $n8 $n1 $n2 $n3 $n4 $n13 $n14 $n15 $n16 $n9 $n10 $n11 $n12 $n17 $n18 $n19 $n20 1 0.0 0.0 $g $rhoi
    } elseif { ($Mat_no) == 2 } {
       element Brick20N $tag $n5 $n6 $n7 $n8 $n1 $n2 $n3 $n4 $n13 $n14 $n15 $n16 $n9 $n10 $n11 $n12 $n17 $n18 $n19 $n20 2 0.0 0.0 $g $rhob 
    } else {
       element Brick20N $tag $n5 $n6 $n7 $n8 $n1 $n2 $n3 $n4 $n13 $n14 $n15 $n16 $n9 $n10 $n11 $n12 $n17 $n18 $n19 $n20 3 0.0 0.0 $g $rhoo
    }
    
    #If this is a pile element, save it in els_in_pile.dat, for replacing after self-weight loading
    if { ($Mat_no) == $PBMat } {
      puts $els_in_pb "$tag "
    }
   
    #Output elements for visualizing
    puts $oelements "$tag $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 $n9 $n10 $n11 $n12 $n13 $n14 $n15 $n16 $n17 $n18 $n19 $n20 $gar $Mat_no"
  
  } else {
    if {($Mat_no) == 1 } {
       element Brick8N $tag $n5 $n6 $n7 $n8 $n1 $n2 $n3 $n4 1 0.0 0.0 $g $rhoi
    } elseif { ($Mat_no) == 2 } {
       element Brick8N $tag $n5 $n6 $n7 $n8 $n1 $n2 $n3 $n4 2 0.0 0.0 $g $rhob
    } else {
       element Brick8N $tag $n5 $n6 $n7 $n8 $n1 $n2 $n3 $n4 3 0.0 0.0 $g $rhoo
    }
    
    #If this is a plastic bowl element, save it in PBElements.dat
    if { ($Mat_no) == $PBMat } {
      puts $els_in_pb "$tag "
    }
    
    #Output elements for visualizing
    puts $oelements "$tag $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 $gar $Mat_no"
  }
    


}
puts "Finished creating $num_of_element elements..."
close $els_in_pb
close $elements
close $oelements


