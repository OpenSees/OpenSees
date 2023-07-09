#================================================================
# Apply_LateralLoad.tcl
# --Create lateral load at pile head
# Zhaohui Yang and Boris Jeremic (@ucdavis.edu)
# March 31, 2002
#================================================================

#Define total lateral load
set Ptot [expr 200.0*$kN];

#Open file nodes_ptop.dat containing nodes on which to apply load
set nd_ptop [open "nodes_ptop.dat" "r"]
gets $nd_ptop num_of_nd_ptop

#Set load per node
set p_pernode [expr $Ptot/$num_of_nd_ptop];

#Apply lateral load at specified top nodes
pattern Plain 2 "Linear" {
   # Read nodes on which lateral load will be applied
   for {set i 0} {$i < $num_of_nd_ptop} {incr i} {
     gets $nd_ptop Toptag
     puts " loading Node $Toptag with Fx = $p_pernode ..."
     #Create the nodal load
     #______tag______Fx______Fy__Fz
     load  $Toptag  $p_pernode 0.0 0.0 -pattern 2
   }

}					      
puts "Finished applying lateral load to specified top nodes..."
close $nd_ptop
