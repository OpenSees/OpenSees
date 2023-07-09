# written: fmk
# date: 02/99
#
# purpose: example1 in g3intro.tex
#
# $Revision: 1.1.1.1 $
# $Date: 2000-09-15 08:23:09 $
# $Source: /usr/local/cvs/OpenSees/EXAMPLES/TclPlaneTruss/example2.tcl,v $


#create the ModelBuilder object
model Tcl2dTruss

# build the model 
node 1   0.0  0.0
node 2 144.0  0.0
node 3 168.0  0.0
node 4  72.0 96.0

material Elastic 1 3000
material ElasticPP 2 3000 .003

truss 1 1 4 10.0 1
truss 2 2 4 5.0 1
truss 3 3 4 5.0 2

# fix the left node, right node on rollar
fix 1 1 1 
fix 2 1 1
fix 3 1 1

pattern Plain 1 Linear {
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
recorder Node Example.out disp -load -nodes 4 -dof 1 2

# perform the analysis
analyze 10

# print the results at node and at all elements
print node 4
print ele
for {set i 1} {$i <= 10} {incr i 1} {playback $i}
