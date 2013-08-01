model BasicBuilder -ndm 2 -ndf 3


set matID [lindex $argv 0]
set motionID [lindex $argv 1]
    
source matProperties.tcl.$matID
set nodes [source model.tcl]
source motion.tcl.$motionID
set node1 [lindex $nodes 0]
set numBay [lindex $nodes 1]
set numFloor [lindex $nodes 2]

set a "recorder EnvelopeNode -file NodeDisp.out.$matID.$motionID -time -node $node1 -dof 1 disp"
eval $a
set a "recorder EnvelopeNode -file NodeAccel.out.$matID.$motionID -time -timeSeries 2 -node $node1 -dof 1 accel"
eval $a

set iNode 1
set jNode [expr $numBay + 2]
set nodesIDrift  ""
set nodesJDrift  ""
set nodeTagsBase ""
for {set i 1} {$i<=$numFloor} {incr i 1} {
    set nodesIDrift "$nodesIDrift $iNode"
    set nodesJDrift "$nodesJDrift $jNode"
    incr iNode [expr $numBay+1]
    incr jNode [expr $numBay+1]
    
}

set a "recorder EnvelopeDrift -file NodeDrift.out.$matID.$motionID -iNode $nodesIDrift -jNode $nodesJDrift -dof 1 -perpDirn 2"
eval $a

for {set i 1} {$i<=[expr $numBay+1]} {incr i 1} {
    set nodeTagsBase "$nodeTagsBase $i"
}

set a "recorder Node -file NodeReaction.out.$matID.$motionID -node $nodeTagsBase -dof 1 reaction"
eval $a


system BandGeneral
constraints Plain
test NormDispIncr 1.0e-12  10 
algorithm Newton
numberer RCM
integrator Newmark  0.5  0.25 
analysis Transient


set t 0.0
set maxT [expr (1+$nPts)*$dT];
set ok 0.0
set maxD 0.0

while {$ok == 0 && $t < $maxT} {
    set ok [analyze 1 $dT]
    if {$ok != 0} {
	test NormDispIncr 1.0e-8  10 1
	algorithm ModifiedNewton -initial
	set ok [analyze 1 $dT]
	test NormDispIncr 1.0e-8  10 
	algorithm Newton
    }
    
    set t [getTime]
}
