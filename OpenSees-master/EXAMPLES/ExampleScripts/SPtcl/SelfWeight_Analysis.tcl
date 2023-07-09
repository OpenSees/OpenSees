#================================================================
# SelfWeight_Analysis.tcl
# --Apply self weight before push-over test
# Zhaohui Yang and Boris Jeremic (@ucdavis.edu)
# March 31, 2002
#================================================================

#Before selfweight analysis, clean previously defined analysis
wipeAnalysis

#Apply Ko consolidation B.C.'s and self-weight
source Apply_KoBC.tcl
source Apply_SelfWeight.tcl

# Define load control parameters
set Steps_sw 10
set lf_sw [expr 1.0/$Steps_sw];

# Define analysis parameters
#set Threshhold 10e-8;
set Threshhold 0.3;
set MaxIte 20;
set Show 1; #parameters indicating show or not show test norm

system UmfPack
constraints Plain
test NormDispIncr $Threshhold $MaxIte $Show
integrator LoadControl $lf_sw 1 $lf_sw $lf_sw
algorithm Newton
numberer RCM
analysis Static

puts "Self weight analysis begins..."
analyze $Steps_sw
puts "Self weight analysis finished..."

