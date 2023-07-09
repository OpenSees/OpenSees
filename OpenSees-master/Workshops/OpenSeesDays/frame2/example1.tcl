
# source in the model
source example.tcl

# create analysis
system BandGeneral
constraints Transformation
numberer RCM
test NormDispIncr 1.0e-12  10 3
algorithm Newton
integrator LoadControl 0.1
analysis Static

#perform analysis
analyze 10

# Print out the state of node 3
print node 3




