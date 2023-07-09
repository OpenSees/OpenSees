# ----------------------------
# Start of model generation
# ----------------------------

# Create ModelBuilder with 3 dimensions and 6 DOF/node
model basic -ndm 3 -ndf 3

# create the material
nDMaterial ElasticIsotropic   1   1000   0.25  6.75 

# Define geometry
# ---------------

# define some  parameters

set mass 0.1

set Quad  quad3d

set eleArgs "PlaneStress2D"
set thick 1.0;

set xWidth 10.0
set yWidth 10.0
set height 5.0

set numX 5; # num elements in x direction
set numY 4; 
set numZ 6;

set nodeNum 1

set incrX [expr $xWidth/(1.0*$numX)]
set incrY [expr $yWidth/(1.0*$numY)]
set incrZ [expr $height/(1.0*$numZ)]

set yLoc 0.0;
set zLoc 0.0;

for {set i 0} {$i <= $numZ} {incr i 1} {
    set xLoc 0.0;
    for {set j 0} {$j <= $numX} {incr j 1} {
	node $nodeNum $xLoc $yLoc $zLoc -mass $mass $mass $mass
	incr nodeNum
	set xLoc [expr $xLoc + $incrX]
    }
    set zLoc [expr $zLoc + $incrZ]
}

set eleTag 1
for {set j 0} {$j <$numZ} {incr j 1} {
    set iNode [expr 1 + $j*($numX+1)]
    set jNode [expr $iNode + 1]
    set kNode [expr $iNode + $numX + 2]
    set lNode [expr $kNode - 1]
    for {set i 0} {$i <$numX} {incr i 1} {
	element quad3d $eleTag $iNode $jNode $kNode $lNode $thick "PlaneStress2D" 1
	incr eleTag; incr iNode; incr jNode; incr kNode; incr lNode;
    }
}

incr nodeNum -1

timeSeries Linear 1
pattern Plain 1 1 {
    load $nodeNum 10.0 0.0 0.0
}

for {set i 1} {$i <= $numX} {incr i 1} {
    fix $i 1 1 1
}

for {set i [expr $numX+1]} {$i <= $nodeNum} {incr i 1} {
    fix $i 0 1 0
}


set xLoc 0.0;
set zLoc 0.0;

incr nodeNum
set nodeYdirnStart $nodeNum

for {set i 0} {$i <= $numZ} {incr i 1} {
    set yLoc 0.0;
    for {set j 0} {$j <= $numY} {incr j 1} {
	node $nodeNum $xLoc $yLoc $zLoc -mass $mass $mass $mass
	incr nodeNum
	set yLoc [expr $yLoc + $incrY]
    }
    set zLoc [expr $zLoc + $incrZ]
}


for {set j 0} {$j <$numZ} {incr j 1} {
    set iNode [expr $nodeYdirnStart + $j*($numY+1)]
    set jNode [expr $iNode + 1]
    set kNode [expr $iNode + $numY + 2]
    set lNode [expr $kNode - 1]
    for {set i 0} {$i <$numY} {incr i 1} {
	element quad3d $eleTag $iNode $jNode $kNode $lNode $thick "PlaneStress2D" 1
	incr eleTag; incr iNode; incr jNode; incr kNode; incr lNode;
    }
}


incr nodeNum -1

for {set i 1; set j $nodeYdirnStart} {$i <= $numY} {incr i 1; incr j 1} {
    fix $j 1 1 1
}

for {set i [expr $nodeYdirnStart+$numY]} {$i <= $nodeNum} {incr i 1} {
    fix $i 1 0 0
}

for {set i [expr 1 + $numX]; set j [expr $nodeYdirnStart+$numY]} {$i < $nodeYdirnStart} {incr i $numX; incr j $numY} {
    remove sp $i 2
    remove sp $j 1
    equalDOF $i $j 1 2 3
}

integrator LoadControl  1.0  1   1.0   10.0
test EnergyIncr  1.0e-12    10         0
algorithm Newton
numberer RCM
constraints Plain 
system ProfileSPD
analysis Static

# Perform the analysis
analyze   10     

# create the display
recorder display g3 10 10 800 200 -wipe
prp -120 -100.0 100.0
vup 0 0 1
viewWindow -30 30 -10 10
display 1 4 5

# --------------------------
# End of recorder generation
# --------------------------

# ---------------------------------------
# Create and Perform the dynamic analysis
# ---------------------------------------

# Remove the static analysis & reset the time to 0.0
wipeAnalysis
setTime 0.0

# Now remove the loads and let the beam vibrate
remove loadPattern 1

# Create the transient analysis
test EnergyIncr  1.0e-12    10         0
algorithm Newton
numberer RCM
constraints Plain 
integrator Newmark 0.5 0.25
system ProfileSPD
#integrator GeneralizedMidpoint 0.50
analysis Transient

# Perform the transient analysis (20 sec)
#       numSteps  dt
analyze  500     0.5

