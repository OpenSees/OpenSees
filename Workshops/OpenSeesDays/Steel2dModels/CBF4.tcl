# EXAMPLE inverted V BRACING

source CBFbase.tcl

# add the inverted V braces
geomTransf Corotational 3
for {set story 1} {$story <= $numStory} {incr story 1} {
    set colLine1 1; set colLine2 2; set colLine3 3;
    set floor1 $story; set floor2 [expr $story +1];
    set theSection [lindex $braceSizes [expr $story -1]]
    # proc HSSbrace {eleTag iNode jNode secType matTag numSeg Im transfTag args}
    HSSbrace $colLine1$floor1$colLine2$floor2 $colLine1$floor1 $colLine2$floor2 $theSection 3 4 $imperfection 3 
    HSSbrace $colLine3$floor1$colLine2$floor2 $colLine3$floor1 $colLine2$floor2 $theSection 3 4 $imperfection 3 
}




