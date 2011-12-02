#create the ModelBuilder object
model basic -ndm 3 -ndf 6

# add nodes - command: node nodeId xCrd yCrd
node 1   0.0  0.0 0.0
node 2   0.0  0.0 0.0

uniaxialMaterial ElasticBilin 1 1000.0 10000.0 0.1
#uniaxialMaterial Elastic 1 1000.0 
element CoupledZeroLength 1 1 2 1 2 1
#element zeroLength 1 1 2 -mat 1 1 -dir 1 2

# set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?
fix 1 1 1 1 1 1 1
fix 2 0 0 1 1 1 1

set P 120.0


# build the components for the analysis object
system BandSPD
#constraints Penalty 1e16 1e16
constraints Plain
integrator LoadControl 0.01
test NormDispIncr 1.0e-12 3 
algorithm Newton
numberer RCM

# create the analysis object 
analysis Static 

# create a Recorder object for the nodal displacements at node 4
recorder Node -file Example.out -load -nodes 4 -dof 1 2 disp

timeSeries Path 1 -dt 1.0 -values {0.0 1.0 0.0 -1.0 0.0}
pattern Plain 1 1 {
    load 2 $P 0.0 0.0 0.0 0.0 0.0
}

analyze 300
remove loadPattern 1

loadConst -time 0.0
print node 2
print ele 1

pattern Plain 1 1 {
    load 2  0.0 $P 0.0 0.0 0.0 0.0
}
test NormDispIncr 1.0e-12 3 0

analyze 300
remove loadPattern 1
loadConst -time 0.0
print node 2
print ele 1

pattern Plain 1 1 {
    set val [expr sqrt($P*$P/2.0)]
    load 2 $val $val 0.0 0.0 0.0 0.0
}

analyze 300
remove loadPattern 1
loadConst -time 0.0
print node 2
print ele 1
