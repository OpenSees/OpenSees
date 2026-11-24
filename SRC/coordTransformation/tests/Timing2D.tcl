# Timing2D.tcl
#
# Computational cost of the 2d transformations Linear, PDelta (default and
# -consistent tangent) and Corotational for two models with the same numbers
# of elements:
#
#   1) FRAME:  a 3-bay moment frame; every column and beam is split into
#              elements of length 1 (storey height 4, bay width 8), and the
#              number of storeys is chosen to give the requested number of
#              elements (40 elements per storey). Gravity load at every floor
#              node plus lateral floor loads increasing with height, scaled so
#              that the Linear roof drift is 5% of the height.
#   2) COLUMN: a single cantilever with ne elements of length 1 (L = ne),
#              loaded at the tip by an end moment M = (1/4) 2 pi EI/L (a
#              quarter circle for the corotational analysis) plus a small
#              transverse force.
#
# Section A = 1, I = 0.01, E = 1e4 for all members. Newton with a banded
# nonsymmetric solver (RCM numbering for the frame), nsteps load steps; each
# case is run nreps times and the fastest time is reported. The output ends
# with a legend and a note on matrix conditioning, which explains why the two
# models need different numbers of iterations (a learning example).
#
# Elements of length 1 are used on purpose: very short, stocky elements
# (length << section radius 0.1) make Newton diverge with a corotational
# transformation.
#
# usage: OpenSees Timing2D.tcl ?neList? ?nsteps? ?nreps?
#   e.g. OpenSees Timing2D.tcl {1000 2000 5000} 10 3

set neList {2000 5000}
set nsteps 10
set nreps  2
if {[llength $argv] > 0} {set neList [lindex $argv 0]}
if {[llength $argv] > 1} {set nsteps [lindex $argv 1]}
if {[llength $argv] > 2} {set nreps  [lindex $argv 2]}

set E 1.0e4
set A 1.0
set I 1.0e-2

# transformations: {label type options}
set transfs {
    {Linear               Linear       ""}
    {PDelta               PDelta       ""}
    {"PDelta -consistent" PDelta       "-consistent"}
    {Corotational         Corotational ""}
}

proc makeElement {tag i j eleType} {
    global E A I
    if {$eleType == "force"} {
        element forceBeamColumn $tag $i $j 1 1
    } else {
        element elasticBeamColumn $tag $i $j $A $E $I 1
    }
}

proc setupAnalysis {nsteps numberer} {
    system BandGeneral
    numberer $numberer
    constraints Plain
    test NormDispIncr 1.0e-6 50 0
    algorithm Newton
    integrator LoadControl [expr 1.0/$nsteps]
    analysis Static
}

# run nsteps load steps; returns {iterations msPerIteration} or {-1 -1}
proc timeAnalysis {nsteps} {
    set its 0
    set t0 [clock microseconds]
    for {set s 0} {$s < $nsteps} {incr s} {
        if {[analyze 1] != 0} {return {-1 -1}}
        incr its [testIter]
    }
    set t1 [clock microseconds]
    return [list $its [expr ($t1-$t0)/1000.0/$its]]
}

# ---------------------------------------------------------------- frame
# nb bays, ns storeys, members split into nsub elements of length 1
# returns {iterations msPerIteration}; lat = lateral load at the roof
proc runFrame {transf opts ns nsteps eleType lat} {
    global E A I
    set nb 3; set nsub 4
    wipe
    model basic -ndm 2 -ndf 3
    set w [expr 2.0*$nsub]
    # column lines: node 1 + c*100000 + j at height j
    for {set c 0} {$c <= $nb} {incr c} {
        for {set j 0} {$j <= [expr $ns*$nsub]} {incr j} {
            node [expr 1+$c*100000+$j] [expr $c*$w] [expr double($j)]
        }
        fix [expr 1+$c*100000] 1 1 1
    }
    eval geomTransf $transf 1 $opts
    if {$eleType == "force"} {
        section Elastic 1 $E $A $I
        beamIntegration Lobatto 1 1 5
    }
    set et 0
    for {set c 0} {$c <= $nb} {incr c} {
        for {set j 0} {$j < [expr $ns*$nsub]} {incr j} {
            makeElement [incr et] [expr 1+$c*100000+$j] [expr 2+$c*100000+$j] $eleType
        }
    }
    set bn 5000000
    for {set s 1} {$s <= $ns} {incr s} {
        set lev [expr $s*$nsub]
        for {set b 0} {$b < $nb} {incr b} {
            set prev [expr 1+$b*100000+$lev]
            for {set k 1} {$k < [expr 2*$nsub]} {incr k} {
                node [incr bn] [expr $b*$w + $k] [expr double($lev)]
                makeElement [incr et] $prev $bn $eleType
                set prev $bn
            }
            makeElement [incr et] $prev [expr 1+($b+1)*100000+$lev] $eleType
        }
    }
    pattern Plain 1 Linear {
        for {set s 1} {$s <= $ns} {incr s} {
            set lev [expr $s*$nsub]
            load [expr 1+$lev] [expr $lat*$s/$ns] 0.0 0.0
            for {set c 0} {$c <= $nb} {incr c} {
                load [expr 1+$c*100000+$lev] 0.0 -0.002 0.0
            }
        }
    }
    setupAnalysis $nsteps RCM
    return [timeAnalysis $nsteps]
}

# lateral load giving a Linear roof drift of 5% of the height (not timed)
proc frameLateralLoad {ns} {
    set nsub 4
    runFrame Linear "" $ns 1 elastic 1.0
    set d [nodeDisp [expr 1+$ns*$nsub] 1]
    return [expr 0.05*$ns*$nsub/$d]
}

# ---------------------------------------------------------------- column
proc runColumn {transf opts ne nsteps eleType} {
    global E A I
    wipe
    model basic -ndm 2 -ndf 3
    set L  [expr double($ne)]
    set pi [expr acos(-1.0)]
    for {set n 1} {$n <= [expr $ne+1]} {incr n} {
        node $n [expr ($n-1)*$L/$ne] 0.0
    }
    fix 1 1 1 1
    eval geomTransf $transf 1 $opts
    if {$eleType == "force"} {
        section Elastic 1 $E $A $I
        beamIntegration Lobatto 1 1 5
    }
    for {set e 1} {$e <= $ne} {incr e} {
        makeElement $e $e [expr $e+1] $eleType
    }
    set tip [expr $ne+1]
    set M [expr 0.25*2*$pi*$E*$I/$L]
    set P [expr 0.001*$E*$I/($L*$L)]
    pattern Plain 1 Linear {
        load $tip 0. $P $M
    }
    setupAnalysis $nsteps Plain
    return [timeAnalysis $nsteps]
}

# ---------------------------------------------------------------- output
# best of nreps: cmd is the run command without the transformation name
proc best {nreps cmd} {
    set bestT -1; set its -1
    for {set r 0} {$r < $nreps} {incr r} {
        lassign [eval $cmd] i t
        if {$i < 0} {return {-1 -1}}
        if {$bestT < 0 || $t < $bestT} {set bestT $t; set its $i}
    }
    return [list $its $bestT]
}

# a cell is 3 lines: total time, time per iteration, iterations per step
proc cell {res} {
    global nsteps
    lassign $res its t
    if {$its < 0} {return [list [format "%-19s" "  diverged"] [format "%-19s" ""] [format "%-19s" ""]]}
    return [list [format "%-19s" [format "%9.3f s" [expr $t*$its/1000.0]]] \
                 [format "%-19s" [format "%9.3f ms" $t]] \
                 [format "%-19s" [format "%7.1f" [expr double($its)/$nsteps]]]]
}

proc header {} {
    global transfs
    set hdr [format "%-9s %-8s %-13s" "elements" "element" ""]
    foreach tr $transfs {append hdr [format " | %-19s" [lindex $tr 0]]}
    puts $hdr
    puts [string repeat "-" [string length $hdr]]
}

# print one model/element row: cells is a list of cell results
proc printRow {ne ele cells} {
    set labels {"total" "per iteration" "its per step"}
    for {set k 0} {$k < 3} {incr k} {
        if {$k == 0} {
            set line [format "%-9d %-8s %-13s" $ne $ele [lindex $labels $k]]
        } else {
            set line [format "%-9s %-8s %-13s" "" "" [lindex $labels $k]]
        }
        foreach c $cells {append line " | " [lindex $c $k]}
        puts $line
    }
    puts ""
}

puts "Timing of 2d transformations ($nsteps load steps, best of $nreps runs)"

puts ""
puts "1) FRAME: 3 bays, members split into elements of length 1, 5% roof drift"
puts ""
header
foreach ne $neList {
    set ns [expr max(1, int(round($ne/40.0)))]
    set lat [frameLateralLoad $ns]
    foreach ele {elastic force} {
        set cells {}
        foreach tr $transfs {
            lassign $tr label type opts
            lappend cells [cell [best $nreps [list runFrame $type $opts $ns $nsteps $ele $lat]]]
        }
        printRow [expr 40*$ns] $ele $cells
    }
}

puts "2) COLUMN: one cantilever of elements of length 1, end moment (quarter circle)"
puts ""
header
foreach ne $neList {
    foreach ele {elastic force} {
        set cells {}
        foreach tr $transfs {
            lassign $tr label type opts
            lappend cells [cell [best $nreps [list runColumn $type $opts $ne $nsteps $ele]]]
        }
        printRow $ne $ele $cells
    }
}

puts "Each cell has three lines:"
puts "  total          wall-clock time of the whole analysis ($nsteps load steps), in s"
puts "  per iteration  total time divided by the number of Newton iterations, in ms;"
puts "                 the cost of one tangent assembly + residual + solve"
puts "  its per step   average number of Newton iterations per load step"
puts "                 (total iterations / $nsteps)"
puts ""
puts "Note on matrix conditioning (why the two models need different numbers of"
puts "iterations for the same number of elements):"
puts "  The condition number of the stiffness matrix is roughly the ratio of its"
puts "  largest to its smallest stiffness. The largest belongs to a single short"
puts "  element (axial EA/l); the smallest to the softest way the whole structure"
puts "  can deform. For the COLUMN that is bending of the whole cantilever, about"
puts "  3EI/L^3 with L = n l, so the condition number grows like n^4: about 3e13"
puts "  for 1000, 5e14 for 2000 and ~1e16 for 5000 elements here (it also grows"
puts "  with element length, like l^2). Double precision has about 16 digits, so the"
puts "  solve keeps only a few correct digits, and Newton needs extra iterations just"
puts "  to remove round-off: even the linear problem (Linear, elastic) takes 4"
puts "  iterations per step instead of 2, and with 20000 elements it does not"
puts "  converge at all."
puts "  The FRAME with the same number of elements is short and stiff in comparison"
puts "  (many members in parallel, many supports), so its matrix is well"
puts "  conditioned: Linear needs 2 iterations per step (one solve plus one check),"
puts "  and the extra iterations of PDelta (gravity axial forces) and Corotational"
puts "  (large displacements) are due only to geometric nonlinearity."
puts "  In the frame most of the time per iteration goes into the equation solve,"
puts "  so the transformations differ little; in the column the solve is cheap"
puts "  and the cost of the transformation itself shows more clearly."
puts "  PDelta: by default its tangent omits the coupling of the axial force with"
puts "  the transverse end displacement (the term (Delta/L) dN/du), so Newton"
puts "  converges only linearly and needs more iterations than Corotational in the"
puts "  frame. With -consistent the term is included: the tangent is exact (but"
puts "  nonsymmetric, so use a nonsymmetric solver such as BandGeneral or UmfPack)"
puts "  and Newton converges quadratically again."
