#================================================================
# Pushover_Analysis.tcl
# --Apply lateral load and perform pushover analysis
# Zhaohui Yang and Boris Jeremic (@ucdavis.edu)
# March 31, 2002
#================================================================

# Clean previously defined analysis method
wipeAnalysis


# Apply Pushover B.C.'s and lateral load for Push-over analysis
#source Apply_PushoverBC.tcl
#source Apply_LateralLoad.tcl


# Define load control parameters
set Steps_po 180
set lf_po [expr 1.0/$Steps_po];

# Define analysis parameters
set Threshhold 0.3;
set MaxIte 20;
set Show 1; #parameters indicating show or not show test norm

#system BandGeneral
#system UmfPack
system SparseGeneral -piv 
constraints Plain
test NormDispIncr $Threshhold $MaxIte $Show
integrator LoadControl $lf_po 1 $lf_po $lf_po
algorithm Newton
numberer RCM
analysis Static

puts "Pushover analysis in $Steps_po steps begins..."
analyze $Steps_po
puts "Pushover analysis finished..."
