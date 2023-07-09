#================================================================
# ForceTopPlateMoveDownHorizontally.tcl
# -- Do what the file name says
# Zhaohui Yang UC Davis
# April 6, 2002
#================================================================


#Open file containing shear nodes on which shear is applied
set shearnds [open "shearnodes.dat" "r"]
gets $shearnds line

set num_shearnds [lindex $line 0];
set MidNode      [lindex $line 1];
puts " Top node = $num_shearnds, Middle top node = $MidNode."

#Read in top nodes to apply multiple-node constraints
for {set i [expr $num_of_node-$num_shearnds+1]} {$i < $num_of_node} {incr i} {
   equalDOF $num_of_node $i 3
   #puts " Added MP constraints of Z dir between $i and $num_of_node..."
}
puts " Added MP constraints of Z dir between all top nodes..."
close $shearnds

# Apply reference load on the top
pattern Plain 3 "Linear" {
   load $MidNode 0.0 0.0 [expr -1*$P]
}

