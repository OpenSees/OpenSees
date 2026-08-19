# Regression model for partitioned loads, mass, floating nodes, and equalDOF.
# Run with: mpiexec -n 4 OpenSeesMP metisPartitionState.tcl

proc csv {values} {
    return [join [lsort -integer $values] ,]
}

model BasicBuilder -ndm 2 -ndf 3
geomTransf Linear 1

# Five mesh nodes and four floating nodes exercise both METIS node ownership
# and nodes that do not occur in the element connectivity.
for {set tag 1} {$tag <= 5} {incr tag} {
    node $tag [expr {double($tag - 1)}] 0.0
}
node 100 2.0 1.0
node 101 5.0 0.0
node 102 5.0 1.0
node 103 6.0 0.0

foreach tag {1 2 3 4 5 100 101 102 103} {
    mass $tag [expr {double($tag)}] 0.0 0.0
}

for {set tag 1} {$tag <= 4} {incr tag} {
    element elasticBeamColumn $tag $tag [expr {$tag + 1}] \
        0.02 2.0e11 8.0e-6 1
}

# A mesh-to-floating pair must follow the mesh node on every rank where the
# constraint is needed. A floating-to-floating pair must have one owner.
equalDOF 3 100 1 2 3
equalDOF 101 102 1 2 3
fix 103 1 1 1

timeSeries Linear 1
pattern Plain 1 1 {
    load 100 10.0 0.0 0.0
    load 103 20.0 0.0 0.0
    eleLoad -ele 1 2 3 4 -type -beamUniform -2.0 0.0
}

partition -objective volume -seed 42

set localNodes [lsort -integer [getNodeTags]]
set masses {}
foreach tag $localNodes {
    lappend masses "$tag:[format %.1f [nodeMass $tag 1]]"
}

puts "METIS_STATE rank=[getPID] ranks=[getNP] nodes=[csv $localNodes] elements=[csv [getEleTags]] eleloads=[csv [getEleLoadTags 1]] nodeloads=[csv [getNodeLoadTags 1]] fixed=[csv [getFixedNodes]] constrained=[csv [getConstrainedNodes]] mp100=[csv [getRetainedNodes 100]] mp102=[csv [getRetainedNodes 102]] masses=[join $masses ,]"
