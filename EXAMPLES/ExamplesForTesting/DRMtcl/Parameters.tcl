#================================================================
# Parameters.tcl
# --Set up some general parameters
# Zhaohui Yang UC Davis
# March 31, 2002
#================================================================

set PI [expr 4*atan(1.0)]
source $Dir/Units.tcl

# gravitational constant in m/sec^2
set g [expr -9.81*$m/pow($sec, 2)]	

# number of nodes per element
set nonper 20;

