#================================================================
# PBowl_Analysis.tcl
# --Apply plastic bowl load and perform transient analysis
# Zhaohui Yang & Boris Jeremic UC Davis
# Nov. 1, 2002
#================================================================

# Apply plastic bowl loading
source $Dir/Apply_PlasticBowlLoading.tcl

# Clean previously defined analysis method
wipeAnalysis

#set previous load constant
#setTime 0.0

# Set up recorders for analysis
source $Dir/SetupPBowlAnalysisRecorder.tcl

# Define load control parameters
# changed to 4 (from 121) so that it does it in reasonable time with insure
set Steps_po 4

# Define analysis parameters
set Threshhold 0.01;
set MaxIte 20;
set Show 1; #parameters indicating show or not show test norm

#system BandGeneral
#system UmfPack
system SparseGeneral -piv 
constraints Penalty 1.0e8 1.0e8

test NormDispIncr $Threshhold $MaxIte $Show
#test NormUnbalance $Threshhold $MaxIte $Show

integrator Newmark 0.505 0.25251 0.0343 0.06 0.06 0.06
algorithm Newton
numberer RCM
analysis Transient

puts "Plastic bowl loading analysis in $Steps_po steps begins..."
#analyze $Steps_po
#for {set ipo 0} {$ipo < $Steps_po} {incr ipo} { 
for {set ipo 0} {$ipo < $Steps_po} {incr ipo} { 
    analyze 1 0.02
}
puts "Plastic bowl loading analysis finished..."
