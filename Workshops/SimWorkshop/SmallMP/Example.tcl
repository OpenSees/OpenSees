set pid [getPID]
set numP  [getNP]
set count 0;
source ReadRecord.tcl

foreach gMotion [glob -nocomplain -directory GM *.AT2] {
    if {[expr $count % $numP] == $pid}  {
	source model.tcl
	source analysis.tcl
	
	set ok [doGravity]
	
	loadConst -time 0.0
	
	if {$ok == 0} {
	    set gMotionName [string range $gMotion 0 end-4 ]
	    
	    set g 384.4
	    ReadRecord ./$gMotionName.AT2 ./$gMotionName.dat dT nPts
	    
	    timeSeries Path 1 -filePath $gMotionName.dat -dt $dT -factor $g	    
	    
	    if {$nPts != 0} {
		
		recorder EnvelopeNode -file $gMotionName.out -node 3 4 -dof 1 2 3 disp
		
		doDynamic $dT $nPts

		file delete $gMotionName.dat


		if {$ok == 0} {
		    puts "$gMotionName OK"
		} else {
		    puts "$gMotionName FAILED"
		}
		
	    } else {
		puts "$gMotion - NO RECORD"
	    }
	}

	wipe
    }
    incr count 1;
}



