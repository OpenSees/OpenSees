#================================================================
# Apply_SelfWeight.tcl
# --Create self weight load for all elements
# Zhaohui Yang and Boris Jeremic (@ucdavis.edu)
# March 31, 2002
#================================================================

#Before apply self-weight load, define a linear load pattern
pattern Plain 1 "Linear" {  

  #Open file elements.dat containing element info
  set elements [open "elements.dat" "r"]
  gets $elements num_of_element
  					      
  for {set i 0} {$i < $num_of_element} {incr i} {
    gets $elements line
    set tag  [lindex $line  0];
    #set n1   [lindex $line  1];
    #set n2   [lindex $line  2];
    #set n3   [lindex $line  3];
    #set n4   [lindex $line  4];
    #set n5   [lindex $line  5];
    #set n6   [lindex $line  6];
    #set n7   [lindex $line  7];
    #set n8   [lindex $line  8];
    #set n9   [lindex $line  9];
    #set n10  [lindex $line 10];
    #set n11  [lindex $line 11];
    #set n12  [lindex $line 12];
    #set n13  [lindex $line 13];
    #set n14  [lindex $line 14];
    #set n15  [lindex $line 15];
    #set n16  [lindex $line 16];
    #set n17  [lindex $line 17];
    #set n18  [lindex $line 18];
    #set n19  [lindex $line 19];
    #set n20  [lindex $line 20];
    #set gar  [lindex $line 21];
    #set Mat_no [lindex $line 22];
  
    #Create elemental load-- SelfWeight
    eleLoad -ele $tag -type -BrickW  
  }
  
}
puts "Finished applying self-weight to all elements..."
close $elements
