
# timeSeries Linear
timeSeries Linear 1

#TCL script: Definitions for hinge
logFile "test.log"

set a_sl 1.0
set cunit 1.0

proc LvOverh { M V h } {
	if {abs($V) < 1e-6} {
		return 4.0
	} else {
		return [expr $M / abs($V) / $h]
	}
}
