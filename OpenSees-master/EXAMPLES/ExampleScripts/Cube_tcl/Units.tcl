#================================================================
# Units.tcl
# --Defining units
#
# Note: In OpenSees soil model implemented by UC Davis 
# uses m, N, and kg as basic unit system. But the input
# data can use any unit, as long as you specify it in 
# your input. For example, set x [expr $x*$in] means x 
# is inch. Internally it will converted to meter.
#
# Zhaohui Yang UC Davis
# April 6, 2002
#================================================================

#Basic units
set m   1.0; # meter for length
set sec 1.0; # second for time
set kg  1.0; # Kilogram for mass


#Other units
# angle
set rad 1.0;
set deg [expr $PI/180.0*$rad];

# length
set cm  0.01;
set in  0.0254;
set ft [expr 12.0*$in];

# mass
set lbs 0.4536;
set kip [expr 1000.0*$lbs];
set ton [expr 1000.0*$kg]

# force
set N   [expr $kg * $m / ($sec * $sec)] ;
set kN [expr 1000.0*$N];

# pressure
set Pa  [expr 1.0*$N/pow($m, 2)];
set kPa [expr 1000.0*$Pa];
set pcf [expr $lbs/pow($ft,3)];	# pcf = #/cubic foot
set ksi [expr $kip/pow($in,2)];
set psi [expr $ksi/1000.];
