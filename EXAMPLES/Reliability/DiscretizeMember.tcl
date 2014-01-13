proc DiscretizeMember {ndI ndJ numEle eleType integration transfTag startNodeTag eleTag} {

    upvar $startNodeTag nodeTag

    if {$numEle <= 1} {
	element $eleType $eleTag $ndI $ndJ $transfTag $integration
	incr eleTag
	return $eleTag
    }

    set Xi [nodeCoord $ndI X]
    set Yi [nodeCoord $ndI Y]

    set Xj [nodeCoord $ndJ X]
    set Yj [nodeCoord $ndJ Y]

    set dX [expr ($Xj-$Xi)/$numEle]
    set dY [expr ($Yj-$Yi)/$numEle]

    set threeD 1
    if {[llength [nodeCoord $ndI]] < 3} {
	set threeD 0
    } else {
	set Zi [nodeCoord $ndI Z]
	set Zj [nodeCoord $ndJ Z]
	set dZ [expr ($Zj-$Zi)/$numEle]
    }

    set nodes(0) $ndI
    set nodes($numEle) $ndJ
    
    # Generate interim nodes
    for {set i 1} {$i < $numEle} {incr i} {
	if {$threeD} {
	    node $nodeTag [expr $Xi+$i*$dX] [expr $Yi+$i*$dY] [expr $Zi+$i*$dZ]
	} else {
	    node $nodeTag [expr $Xi+$i*$dX] [expr $Yi+$i*$dY]
	}
	set nodes($i) $nodeTag
	incr nodeTag
    }

    # Generate end elements
    element $eleType $eleTag $ndI $nodes(1) $transfTag $integration
    incr eleTag

    element $eleType $eleTag $nodes([expr $numEle-1]) $ndJ $transfTag $integration
    incr eleTag

    # Generate interior elements
    for {set i 1} {$i < [expr $numEle-1]} {incr i} {
	element $eleType $eleTag $nodes($i) $nodes([expr $i+1]) $transfTag $integration
	incr eleTag
    }

    return $eleTag
}
