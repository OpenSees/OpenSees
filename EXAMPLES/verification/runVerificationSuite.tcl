# script to run all the verification scripts
# PASS/FAILURE results in file results.out when run

# open results file such that it is cleared out of any data
set results [open results.out w]
close $results

source PlanarTruss.tcl
source PortalFrame2d.tcl
source EigenFrame.tcl
source EigenFrame.Extra.tcl

