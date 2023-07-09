# ReadRecord.tcl
# ------------------------------------------------------------------------------------------------------------
#
# Written: fmk
# Date: July 2010

# A procedure which parses a ground motion record from the PEER
# strong motion database by finding dt in the record header, then
# echoing data values to the output file.
#
# Formal arguments
#	inFilename -- file which contains PEER strong motion record
#	outFilename -- file to be written in format G3 can read
#	dt -- time step determined from file header
#	nPts -- number of data points from file header
#
# Assumptions
#	The header in the PEER record is, e.g., formatted as 1 of following:
#  1) new PGA database
#	 PACIFIC ENGINEERING AND ANALYSIS STRONG-MOTION DATA
#	  IMPERIAL VALLEY 10/15/79 2319, EL CENTRO ARRAY 6, 230                           
#	  ACCELERATION TIME HISTORY IN UNITS OF G                                         
#	  3930 0.00500 NPTS, DT

#   2) old SMD database
#	 PACIFIC ENGINEERING AND ANALYSIS STRONG-MOTION DATA
#	  IMPERIAL VALLEY 10/15/79 2319, EL CENTRO ARRAY 6, 230                           
#	  ACCELERATION TIME HISTORY IN UNITS OF G                                         
#	  NPTS=  3930, DT= .00500 SEC


proc ReadRecord {inFilename outFilename dt nPts} {

    # Pass dt by reference
    upvar $dt DT
    upvar $nPts NPTS
    
    # Open the input file and catch the error if it can't be read
    if [catch {open $inFilename r} inFileID] {
	puts stderr "Cannot open $inFilename for reading"
    } else {
	# Open output file for writing
	set outFileID [open $outFilename w]
	
	# Flag indicating dt is found and that ground motion
	# values should be read -- ASSUMES dt is on last line
	# of header!!!
	set flag 0
	
	# Look at each line in the file
	foreach line [split [read $inFileID] \n] {
	    
	    if {[llength $line] == 0} {
		# Blank line --> do nothing
		continue
	    } elseif {$flag == 1} {
		# Echo ground motion values to output file
		puts $outFileID $line
	    } else {
	
		# Search header lines for dt
		set lengthLine [llength $line]

		if {$lengthLine >= 4} {
		
		    set word0 [lindex $line 0]
		    set wordN [lindex $line [expr $lengthLine-1]]
		    
		    if {$word0 == "NPTS=" } {
			# old SMD format
			set count 0
			foreach word [split $line] {
			    incr count 1
			    if {$word != ""} {
				# Read in the time step
				if {$flag == 1} {
				    set DT $word
				    break;
				}
				if {$flag == 2} {
				    set NPTS [string trim $word ","]
				    set flag 0
				}
				# Find the desired token and set the flag
				if {[string match $word "DT="] == 1} {
				    set flag 1
				}
				# Find the desired token and set the flag
				if {[string match $word "NPTS="] == 1} {
				    set flag 2
				}
			    }
			}

		    } elseif {$wordN == "DT"} {

			# new NGA format
			set count 0;
			foreach word [split $line] {
			    if {$word != ""} {
				if {$count == 0} {
				    set NPTS $word;
				} elseif {$count == 1} {
				    set DT $word;
				} elseif {[string match $word "DT"] == 1} {
				    set flag 1;
				    break;
				}
				incr count 1
			    }
			}
		    }
		}
	    }
	}
	# Close the output file
	close $outFileID
	
	# Close the input file
	close $inFileID
    }
}

