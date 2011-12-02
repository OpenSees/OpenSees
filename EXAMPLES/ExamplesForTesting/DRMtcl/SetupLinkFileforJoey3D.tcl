#================================================================
# SetupLinkFileforJoey3D.tcl
# --Set up link file for visualizer Joey3D
# Zhaohui Yang UC Davis
# March 31, 2002
#================================================================

#File to link output with Joey3D
set head [open "headers.dat" "w"];

#Write number of steps in Pushover analysis
puts $head "$isw $ipo"; 
puts $head "GaussPoint";
puts $head "Displ";
puts $head "Plast";
puts $head "Sigma";
puts $head ".out";


