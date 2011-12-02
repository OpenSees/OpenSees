#================================================================
# Units.tcl
# --Defining units
# Zhaohui Yang UC Davis
# March 31, 2002
#================================================================

#Basic units
set m   1.0; 
set sec 1.0; 
set kg  1.0; 


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

# pcf = #/cubic foot
set pcf [expr $lbs*9.81/pow($ft,3)];	

set ksi [expr $kip*9.81/pow($in,2)];
set psi [expr $ksi/1000.];
