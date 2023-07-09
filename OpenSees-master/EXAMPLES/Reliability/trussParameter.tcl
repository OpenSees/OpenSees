model basic -ndm 2 -ndf 2
#reliability

set L 20.0
set E 30000.0
set A 5.0
set P 25.0

node 1 0.0 0.0; fix 1 1 1
node 2 [expr 0.5*$L] 0.0; fix 2 0 1
node 3 $L 0.0; fix 3 0 1

uniaxialMaterial Elastic 1 $E

element truss 1 1 2 $A 1
element truss 2 2 3 $A 1

pattern Plain 1 Linear {
    load 3 $P 0.0
}

# Material property of an element
parameter 1 element 1 E Area One
addToParameter 1 element 2 E 

# Element property
parameter 2 element 2 A
addToParameter 2 element 1 A

# Nodal load
parameter 3 loadPattern 1 loadAtNode 3 1

# Nodal coordinate
parameter 4 node 3 coord 1

analysis Static

analyze 1

set U [nodeDisp 3 1]
puts "Truss deflection before update: $U"


# Perturb parameters relative to their original values
updateParameter 1 [expr 0.95*$E]
updateParameter 2 [expr 1.05*$A]
updateParameter 3 [expr 1.05*$P]
updateParameter 4 [expr 0.95*$L]

analyze 1

set U [nodeDisp 3 1]
puts "Truss deflection after update: $U"
