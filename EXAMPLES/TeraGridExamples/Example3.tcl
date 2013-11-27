
set pid [getPID]
set numP  [getNP]
set recordsFileID [open "peerRecords.txt" r]
set count 0;

foreach gMotion [split [read $recordsFileID] \n] {
    if {[expr $count % $numP] == $pid}  {

	source model.tcl
	source analysis.tcl
	
	set ok [doGravity]

	loadConst -time 0.0
	
	if {$ok == 0} {
	    set gMotionList [split $gMotion "/"]
	    set gMotionDir  [lindex $gMotionList end-1]
	    set gMotionNameInclAT2 [lindex $gMotionList end]
	    set gMotionName [string range $gMotionNameInclAT2 0 end-4 ]
	    
	    set Gaccel "PeerDatabase $gMotionDir $gMotionName -accel 384.4 -dT dT -nPts nPts"
	    pattern UniformExcitation 2 1 -accel $Gaccel
	    
	    if {$nPts != 0} {
		
		recorder EnvelopeNode -file $gMotionDir$gMotionName.out -node 3 4 -dof 1 2 3 disp
		
		doDynamic $dT $nPts
		
		if {$ok == 0} {
		    puts "$gMotionDir $gMotionName OK"
		} else {
		    puts "$gMotionDir $gMotionName FAILED"
		}
		
	    } else {
		puts "$gMotion - NO RECORD"
	    }
	}

	wipe
    }
    incr count 1;
}



