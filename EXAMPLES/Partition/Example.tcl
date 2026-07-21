# OpenSeesMP partition example: serial, METIS, custom, and samePart.
#
# A small 3D solid pier is solved four ways. Each run prints element ownership,
# the first three eigenvalues, and top-node displacement after gravity and after
# a short transient pulse. The numerical results should agree independent of
# how the mesh is partitioned.
#
# Run from this directory (any process count from 2 to 8 works for the
# parallel modes):
#   OpenSees Example.tcl serial
#   mpiexec -n 4 OpenSeesMP Example.tcl metis
#   mpiexec -n 4 OpenSeesMP Example.tcl custom
#   mpiexec -n 4 OpenSeesMP Example.tcl samePart
#
# Optional smaller mesh for a quick run:
#   OpenSees Example.tcl serial 4 4 20
#   mpiexec -n 4 OpenSeesMP Example.tcl metis 4 4 20

wipe

set pid [getPID]
set np  [getNP]
set mode [expr {$argc > 0 ? [string tolower [lindex $argv 0]] : "serial"}]

set nx 10
set ny 10
set nz 100
if {$argc == 4} {
    foreach {nx ny nz} [lrange $argv 1 3] break
    foreach n [list $nx $ny $nz] {
        if {![string is integer -strict $n] || $n < 1} {
            puts "ERROR: nx, ny, and nz must be positive integers"
            exit 1
        }
    }
} elseif {$argc != 1} {
    puts "usage: Example.tcl serial|metis|custom|samePart ?nx ny nz?"
    exit 1
}

if {$mode ni {serial metis custom samepart}} {
    puts "usage: Example.tcl serial|metis|custom|samePart ?nx ny nz?"
    exit 1
}
if {$mode eq "serial" && $np != 1} {
    puts "ERROR: serial mode requires one process"
    exit 1
}
if {$mode ne "serial" && $np < 2} {
    puts "ERROR: parallel modes in this example require at least 2 MPI processes"
    exit 1
}

# Same concrete-like 2 m x 2 m x 10 m pier as EXAMPLES/LargeMP.
set nn [expr {($nx + 1)*($ny + 1)*($nz + 1)}]
set ne [expr {$nx*$ny*$nz}]
set topNode $nn

model BasicBuilder -ndm 3 -ndf 3
set rho 2.4
nDMaterial ElasticIsotropic 1 2.5e7 0.20 $rho

block3D $nx $ny $nz 1 1 stdBrick 1 {
    1  -1.0  -1.0   0.0
    2   1.0  -1.0   0.0
    3   1.0   1.0   0.0
    4  -1.0   1.0   0.0
    5  -1.0  -1.0  10.0
    6   1.0  -1.0  10.0
    7   1.0   1.0  10.0
    8  -1.0   1.0  10.0
}

# Element mass comes from material rho. Add 50% extra nonstructural mass,
# distributed uniformly over all nodes, to exercise nodal-mass ownership too.
set structuralMass [expr {$rho*2.0*2.0*10.0}]
set nodalMass [expr {0.5*$structuralMass/double($nn)}]
for {set i 1} {$i <= $nn} {incr i} {
    mass $i $nodalMass $nodalMass $nodalMass
}
fixZ 0.0 1 1 1

# Mesh-independent gravity from all element + nodal mass. UniformExcitation
# applies F = -M*ugddot, so positive g acts downward in global Z (direction 3).
timeSeries Constant 1
pattern UniformExcitation 1 3 -accel 1 -fact 9.81

if {$mode eq "metis"} {
    partition
} elseif {$mode eq "custom"} {
    # Explicit element -> part map: contiguous, balanced blocks over all ranks.
    set customPairs {}
    for {set e 1} {$e <= $ne} {incr e} {
        lappend customPairs $e [expr {($e - 1)*$np/$ne}]
    }
    partition -customPartition $ne {*}$customPairs
} elseif {$mode eq "samepart"} {
    # METIS still partitions, but first/last (nonadjacent) elements stay together.
    partition -samePart 2 1 $ne
}

set localEles [lsort -integer [getEleTags]]
set eleSample [lrange $localEles 0 11]
puts "PARTITION mode=$mode rank=$pid count=[llength $localEles] sample=$eleSample"

if {$mode eq "samepart"} {
    set hasFirst [expr {[lsearch -exact $localEles 1] >= 0}]
    set hasLast [expr {[lsearch -exact $localEles $ne] >= 0}]
    if {$hasFirst != $hasLast} {
        puts "ERROR: -samePart failed: elements 1 and $ne are split on rank $pid"
        exit 1
    }
}

if {$np > 1} {
    set numbererType ParallelPlain
    set systemType Mumps
} else {
    set numbererType Plain
    set systemType UmfPack
}

constraints Plain
numberer $numbererType
system $systemType
test NormUnbalance 1.0e-10 20
algorithm Newton
integrator LoadControl 0.1
analysis Static

set ok [analyze 10]
if {$ok != 0} {
    puts "ERROR: gravity analysis failed on rank $pid, code=$ok"
    exit 1
}

# Only the rank retaining the top node prints its displacement.
set ownsTop [expr {[lsearch -exact [getNodeTags] $topNode] >= 0}]
if {$ownsTop} {
    puts [format "RESULT mode=%s stage=gravity node=%d disp={%.16e %.16e %.16e}" \
        $mode $topNode [nodeDisp $topNode 1] [nodeDisp $topNode 2] \
        [nodeDisp $topNode 3]]
}

# All MPI ranks participate in eigen and receive the same eigenvalues.
set lambdas [eigen 3]
if {$pid == 0} {
    puts "RESULT mode=$mode stage=eigen lambda=$lambdas"
}

# Hold gravity constant and apply a one-step horizontal pulse at the top.
set period1 [expr {2.0*acos(-1.0)/sqrt([lindex $lambdas 0])}]
set dt [expr {$period1/20.0}]
loadConst -time 0.0
timeSeries Path 2 -dt $dt -values {0.0 1.0 0.0} -factor 100.0
pattern Plain 2 2 {
    if {$ownsTop} {
        load $topNode 1.0 0.0 0.0
    }
}

wipeAnalysis
constraints Plain
numberer $numbererType
system $systemType
test EnergyIncr 1.0e-10 20
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient

set ok [analyze 40 $dt]
if {$ok != 0} {
    puts "ERROR: transient analysis failed on rank $pid, code=$ok"
    exit 1
}
if {$ownsTop} {
    puts [format "RESULT mode=%s stage=transient node=%d disp={%.16e %.16e %.16e}" \
        $mode $topNode [nodeDisp $topNode 1] [nodeDisp $topNode 2] \
        [nodeDisp $topNode 3]]
}

wipe
exit 0
