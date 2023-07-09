# written: fmk
# date: 02/99
#
# purpose: example2 in g3intro.tex modified for TclModelBuilder
#  and used to show ZeroLength element at node 0
#
# $Revision: 1.2 $
# $Date: 2001-08-15 21:20:53 $
# $Source: /usr/local/cvs/OpenSees/EXAMPLES/ExamplesForTesting/Test.example3.ops,v $


#create the ModelBuilder object
model BasicBuilder -ndm 2 -ndf 2

# build the model 
node 1   0.0  0.0
node 2 144.0  0.0
node 3 168.0  0.0
node 4  72.0 96.0
node 5   0.0  0.0

uniaxialMaterial Elastic 1 3000
uniaxialMaterial ElasticPP 2 3000 .003

element truss 1 1 4 10.0 1
element truss 2 2 4 5.0 1
element truss 3 3 4 5.0 2

uniaxialMaterial Elastic 3 3.0e8
uniaxialMaterial Elastic 4 3.0e2

element zeroLength 4 1 5 -mat 4 3 -dir 1 2 -orient 1 0 0 0 1 0 

#element zeroLength 4 1 5 -mat 3 3 -dirn 1 2 
#element zeroLength  4 1 5 -mat 3 -dirn 1
#element zeroLength  5 1 5 -mat 3 -dirn 2

# fix the left node, right node on rollar
fix 5 1 1 
fix 2 1 1
fix 3 1 1

pattern Plain 1 "Linear" {
    load 4 100 -50
}

# build the components for the analysis object
system BandSPD
constraints Plain
integrator LoadControl 0.1 1 0.1 0.1
test NormUnbalance 1.0e-6
algorithm Newton
numberer RCM

# create the analysis object 
analysis Static 

# create a Recorder object for the nodal displacements at node 4
recorder Node Node.out disp -load -nodes 4 -dof 1 2

# perform the 10 analysis steps
analyze 10

# print the results at node and at all elements
print node 4
print ele
for {set i 1} {$i <= 10} {incr i 1} {playback $i}

wipe
