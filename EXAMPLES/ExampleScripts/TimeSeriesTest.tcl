
model basic -ndm 2 -ndf 2

set PI [expr 2*asin(1.0)]
set PIdiv2 [expr $PI/2.0]

set filePath [open filePath w]
set timePath [open fileTime w]
set both [open fileBoth w]
for {set i 0} {$i <= 11} {incr i 1} {
    puts $filePath $i
    puts $timePath [expr $i/2.0]
    puts $both "$i [expr $i/2.0]"
}

close $filePath
close $timePath
close $both

set records [searchPeerNGA -fault San -magLo 6.5 -magHi 6.75]
    
set seriesCommands [list]
lappend seriesCommands "Linear 1 -factor 2.0"
lappend seriesCommands "Constant 1 -factor 2.0"
lappend seriesCommands "Rectangular 1 2.0 9.0"
lappend seriesCommands "Rectangular 1 2.0 9.0 -factor 2.0"
lappend seriesCommands "Trig 1 1.0 9.0 8"
lappend seriesCommands "Trig 1 1.0 9.0 8 -factor 2.0"
lappend seriesCommands "Trig 1 1.0 9.0 8 -factor 2.0 -shift $PIdiv2"
lappend seriesCommands "Trig 1 1.0 9.0 8 -factor 2.0"
lappend seriesCommands "Trig 1 1.0 9.0 8 -factor 2.0 -shift $PIdiv2"
lappend seriesCommands "Path 1 -dt 1.0 -filePath filePath"
lappend seriesCommands "Path 1 -dt 1.0 -filePath filePath -factor 2.0"
lappend seriesCommands "Path 1 -fileTime fileTime -filePath filePath -factor 2.0"
lappend seriesCommands "Path 1 -dt 1.0 -values {0.0 1.0 2.0 3.0} -factor 2.0"
lappend seriesCommands "PeerNGAMotion 1 [lindex $records 0] 1.0 -dT dt -NPTS nPts"

#foreach record $records {
#    lappend seriesCommands "PeerNGAMotion $record 1.0 -dT dt -NPTS nPts"
#}

foreach seriesCommand $seriesCommands {
    
    # model                                                                         
    set periodStruct 1.0;
    set L 1.0
    set A 1.0;
    set m 1.0;
    
    #derived quantities                                                             
    set PI 3.14159

    set fStruct [expr 1.0/$periodStruct]
    set wStruct [expr 2.0 * $PI / $periodStruct]
    set K [expr $wStruct * $wStruct * $m]
    set E [expr $L * $K / $A]
    
    wipe
    model basic -ndm 1 -ndf 1
    
    node  1       0.0
    node  2       $L -mass $m
    uniaxialMaterial Elastic 1 $E
    element truss 1 1 2 $A 1
    
    fix 1 1

    eval "timeSeries $seriesCommand"

    pattern Plain 1 1 {
	load  2   1.0  
    }

    integrator LoadControl 1.0
    constraints Plain
    system BandGen
    algorithm Linear
    numberer Plain
    constraints Plain
    analysis Static

    recorder Node -file D.out -timeSeries 1 -time -node 2 -dof 1 disp
    recorder Node -file V.out -timeSeries 1 -time -node 2 -dof 1 vel
    recorder Node -file A.out -timeSeries 1 -time -node 2 -dof 1 accel

    puts "$seriesCommand"    
    for {set i 1} {$i <= 10} {incr i 1} {
	analyze 1
	set node2Disp [nodeDisp 2 1]	
	puts "$i $node2Disp"
    }
}

