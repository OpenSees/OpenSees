foreach record {el_centro IELC180 ALH360} {

    wipe

    # set some variables
    set Tn 1.0
    set K 4.0
    set dampRatio 0.02

    #set some constants                                                             
    set g 386.4
    set PI [expr 2.0 * asin(1.0)]

    #derived quantaties
    set Wn [expr 2.0 * $PI / $Tn]
    set M [expr $K / ($Wn * $Wn)]
    set c [expr 2.0*$M*$Wn*$dampRatio]

    # create the model
    wipe
    model basic -ndm 1 -ndf 1

    node  1       0.0
    node  2       1.0  -mass $M

    fix 1 1

    uniaxialMaterial Elastic  1 $K 0.0
    uniaxialMaterial Elastic  2 0.0 $c
    uniaxialMaterial  Parallel 3 1 2

    element truss 1 1 2 1.0 3

    set dT 0.0;
    set nPts 0;
    set factor $g

    source ReadRecord.tcl
    ReadRecord $record.AT2  $record.dat dT nPts	    

    timeSeries Path 1 -filePath $record.dat -dt $dT -factor $g
    pattern UniformExcitation  1 1 -accel 1

    # create the analysis
    constraints Plain
    integrator Newmark 0.5 [expr 1.0/6.0]
    system ProfileSPD
    test NormUnbalance 1.0e-12  6 0
    algorithm Newton
    numberer RCM
    analysis Transient

    set t 0.0
    set maxT [expr (1+$nPts)*$dT];
    set ok 0.0
    set maxD 0.0

    recorder Node -file accel$record.out -timeSeries 1 -node 2 -dof 1 accel
    recorder EnvelopeNode -file disp$record.out -timeSeries 1 -node 2 -dof 1 disp


    while {$ok == 0 && $t < $maxT} {
	set ok [analyze 1 $dT]
	if {$ok != 0} {
	    test NormDispIncr 1.0e-12  100 0
	    algorithm ModifiedNewton -initial
	    set ok [analyze 1 .01]
	    test NormDispIncr 1.0e-12  10 
	    algorithm Newton
	}
    
	set time [getTime]

	set t [expr $t + $dT]
    }
}
