# script to run all the verification scripts
# PASS/FAILURE results in file results.out when run

# open results file such that it is cleared out of any data
set results [open results.out w]
close $results

source sdofTransient.tcl
source PlanarTruss.tcl
source PlanarTruss.Extra.tcl
source PortalFrame2d.tcl
source EigenFrame.tcl
source EigenFrame.Extra.tcl
source AISC25.tcl
source PlanarShearWall.tcl
source PinchedCylinder.tcl

exit
