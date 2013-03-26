set pid [getPID]
set numP  [getNP]
set count 0;
source ReadRecord.tcl

set tStart [clock clicks -milliseconds]

foreach gMotion [glob -nocomplain -directory GM *.AT2] {
    if {[expr $count % $numP] == $pid}  {

	puts "$pid $count $gMotion"

	source model.tcl
	source analysis.tcl
	
	set ok [doGravity]
	
	loadConst -time 0.0
	
	if {$ok == 0} {
	    set gMotionName [string range $gMotion 0 end-4 ]
	    
	    set g 384.4
	    ReadRecord ./$gMotionName.AT2 ./$gMotionName.dat dT nPts
	    
	    timeSeries Path 1 -filePath ./$gMotionName.dat -dt $dT -factor $g	    
	    pattern UniformExcitation 2 1 -accel 1

	    if {$nPts != 0} {
		
		recorder Node -file $gMotionName.out -node 3 4 -dof 1 2 3 disp
		
		doDynamic $dT $nPts
#		doDynamic $dT 2

#		file delete $gMotionName.dat


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


set tEnd [clock clicks -milliseconds]
set duration [expr $tEnd-$tStart]
if {$pid == 0} {
    puts "Duration $duration"
}