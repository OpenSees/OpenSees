# source in the model
source example1.tcl

#set loads constant & reset time in domain to be 0.0
loadConst -time -0.01

# set some parameters
set g 386.4

set motion IELC180
# Source in TCL proc to read PEER SMD record
source ReadSMDFile.tcl

# Read the file
ReadSMDFile $motion.AT2 $motion.acc dt
ReadSMDFile $motion.DT2 $motion.disp dt

# create the load pattern
set accelSeries "Path -filePath $motion.acc -dt $dt -factor $g"
pattern UniformExcitation  2   1  -accel $accelSeries

# set the rayleigh damping factors for nodes & elements
rayleigh 0.0 0.0 0.0 0.0

wipeAnalysis

system ProfileSPD
constraints Plain
test NormUnbalance 1.0e-14  10 
algorithm Newton
numberer RCM
integrator Newmark  0.5  0.25 

analysis Transient

recorder Node -time -file example3.out -node 3 -dof 1 disp

puts "Eigenvalue: [eigen 1]"

analyze 4000 $dt

puts "Eigenvalue: [eigen 1]"

puts [eigen 1]
print node 3