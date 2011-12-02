#================================================================
# IsotropicCompression_Analysis.tcl
# --Isotropic compression analysis before shear analysis
# Zhaohui Yang UC Davis
# April 6, 2002
#================================================================

#Before isotropic compression analysis, clean previously defined analysis
wipeAnalysis

#Apply isotropic pressure
source Apply_IsotropicPressure.tcl

# Set up recorders for Self Weight analysis
source SetupIsotropicCompressionRecorder.tcl 

# Define load control parameters
set Steps_iso 10
set lf_iso [expr 1.0/$Steps_iso];

# Define analysis parameters
#set Threshhold 10e-8;
set Threshhold 1.0e-08;
set MaxIte 10;
set Show 1; #parameters indicating to show or not show residual value

system SparseGeneral
#system UmfPack
constraints Plain
test NormDispIncr $Threshhold $MaxIte $Show
integrator LoadControl $lf_iso 1 $lf_iso $lf_iso
algorithm Newton
numberer RCM
analysis Static

puts "Isotropic consolidation analysis begins..."
#analyze $Steps_sw
for {set iso 0} {$iso < $Steps_iso} {incr iso} { 
    analyze 1
}
puts "Isotropic consolidation analysis finished..."

