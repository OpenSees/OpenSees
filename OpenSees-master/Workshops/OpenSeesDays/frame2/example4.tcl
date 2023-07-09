# source in the model
source example1.tcl

#set loads constant & reset time in domain to be 0.0
loadConst -time -0.01

# set some parameters
set motion IELC180

# Source in TCL proc to read PEER SMD record
source ReadSMDFile.tcl

# Read the file
ReadSMDFile $motion.DT2 $motion.disp dt

# create the load pattern
set dispSeries "Path -filePath $motion.disp -dt $dt -factor 0.3937"

pattern MultiSupport  2   {
    groundMotion 5 Plain -disp $dispSeries
    imposedMotion 1 1 5
    imposedMotion 2 1 5
}

# set the rayleigh damping factors for nodes & elements
rayleigh 0.0 0.0 0.0 0.0

remove sp 1 1
remove sp 2 1

#create the analysis, start by getting rid of old static
wipeAnalysis
system ProfileSPD
constraints Transformation
test NormDispIncr 1.0e-14  10 
algorithm Newton
numberer RCM
integrator Newmark  0.5  0.25 

analysis Transient

recorder Node -time -file example4.out -node 3 -dof 1 disp

puts "Eigen before Gravity: [eigen 1]"

analyze 4000 $dt

puts "Eigen after Gravity: [eigen 1]"
print node 3


