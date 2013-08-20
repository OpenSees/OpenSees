# model builder
model Basic -ndm 2 -ndf 2

# some problem parameters
set L 40.0
set H 10.0
set thick 2.0

set P 10
set nX 9; # numNodes x dirn
set nY 3; # numNodes y dirn

# create material
nDMaterial ElasticIsotropic 1 1000 0.25  3.0

# create nodes
set nodeTag 1
set yLoc 0.0;
for {set i 0} {$i < $nY} {incr i 1} {
    set xLoc 0.0;
    for {set j 0} {$j < $nX} {incr j 1} {
	node $nodeTag $xLoc $yLoc
        set xLoc [expr $xLoc+ $L/($nX-1.0)]
	incr nodeTag
    }
    set yLoc [expr $yLoc+ $H/($nY-1.0)]
}

# boundary conditions
fix  1  1 1
fix $nX 1 1

# create elements
set eleTag 1
for {set i 1} {$i < $nY} {incr i 1} {
    set iNode [expr 1+($i-1)*$nX]; 
    set jNode [expr $iNode+1];
    set kNode [expr $jNode+$nX]
    set lNode [expr $iNode+$nX]
    for {set j 1} {$j < $nX} {incr j 1} {
	element quad $eleTag $iNode $jNode $kNode $lNode $thick "PlaneStress" 1
	incr eleTag
	incr iNode; incr jNode; incr kNode; incr lNode
    }
}

# apply loads
set midNode [expr ($nX+1)/2]
timeSeries Linear 1
pattern Plain 1 1 {
    load $midNode 0 -$P
    load [expr $midNode + $nX*($nY-1)] 0 -$P
}

# create analysis
numberer RCM
algorithm Linear
integrator LoadControl 1
constraints Plain
system SparseSYM
analysis Static

# analyze
analyze 1

print node $midNode

recorder display "Simply Supported Beam" 10     10      800     200    -wipe  
prp 20 5.0 1.0;                                      # projection reference point (prp); defines the center of projection (viewer eye)
vup  0  1 0;                                         # view-up vector (vup) 
vpn  0  0 1;                                         # view-plane normal (vpn)     
viewWindow -30 30 -10 10;                            # coordiantes of the window relative to prp  
display 10 0 5;                                      # the 1st arg. is the tag for display mode
