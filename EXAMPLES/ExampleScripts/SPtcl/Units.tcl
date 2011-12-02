#================================================================
# Units.tcl
# --Defining units
# Zhaohui Yang and Boris Jeremic (@ucdavis.edu)
# March 31, 2002
#================================================================

#Basic units
set m   1.0; # meter for length
set sec 1.0; # second for time
set kg  1.0; # Kilogram for mass
set N   1.0; # Newton for force

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
set kN [expr 1000.0*$N];

# pressure
set Pa  [expr 1.0*$N/pow($m, 2)];
set kPa [expr 1000.0*$Pa];
set pcf [expr $lbs/pow($ft,3)];	# pcf = #/cubic foot
set ksi [expr $kip/pow($in,2)];
set psi [expr $ksi/1000.];
