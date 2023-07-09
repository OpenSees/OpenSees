# first source in the model
source RCFrame.tcl

# Create the gravity loads
set W 4000.0;
timeSeries Linear 1
pattern Plain 1 1 {
    eleLoad -ele 3 -type -beamUniform [expr -$W/$width]
}

# create the analysis
system BandGeneral
constraints Transformation
numberer RCM
test NormDispIncr 1.0e-12  10 3
algorithm Newton
integrator LoadControl 0.1
analysis Static

# perform the analysis
analyze 10

# Print out the state of node 3
print node 3

# print state of element 1
print ele 1

