#================================================================
# Pushover_Analysis.tcl
# --Apply lateral load and perform pushover analysis
# Zhaohui Yang UC Davis
# April 6, 2002
#================================================================

# Clean previously defined analysis method
wipeAnalysis

#set previous load constant
setTime 0.0
loadConst

# Apply uniform axial loading B.C.'s and set up recorders
source ForceTopPlatenMoveDownHorizontally.tcl
source SetupAxialLoadRecorder.tcl

# Define displacement control parameters
set TotDis -0.05
set Steps_al 100
set ndz [expr $TotDis/$Steps_al]

# Define analysis parameters
set Threshhold 1.0e-04;
set MaxIte 20;
set Show 1; #parameters indicating show or not show test norm
set P1 1e12; #Penalty constant
set P2 1e12; #Penalty constant

#system BandGeneral
#system UmfPack
system SparseGeneral -piv 
#constraints Penalty $P1 $P2
constraints Plain
test NormDispIncr $Threshhold $MaxIte $Show
integrator DisplacementControl $MidNode 3 $ndz 1 $ndz $ndz
#integrator DisplacementControl $MidNode 3 $ndz 10 [expr $ndz/10] [expr $ndz*10]
algorithm Newton
numberer RCM
analysis Static

puts "Axial loading analysis in $Steps_al steps begins..."
#analyze $Steps_po
#for {set ipo 0} {$ipo < $Steps_po} {incr ipo} { 
for {set ial 0} {$ial < $Steps_al} {incr ial} { 
    analyze 1
}
puts "Axial loading analysis finished..."
