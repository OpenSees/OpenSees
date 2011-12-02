#================================================================
# Parameters.tcl
# --Set up some general parameters
# Zhaohui Yang
# March 31, 2002
#================================================================

set PI [expr 4*atan(1.0)];
source Units.tcl
set g [expr -9.81*$m/pow($sec, 2)];	# gravitational constant in m/sec^2
set nonper 8; # number of nodes per element
