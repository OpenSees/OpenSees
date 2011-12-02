# written: fmk
# date: 02/99
#
# purpose: example1 in OpenSeesIntro.tex
#
# $Revision: 1.2 $
# $Date: 2000-12-16 06:24:06 $
# $Source: /usr/local/cvs/OpenSees/EXAMPLES/ExampleScripts/example1.tcl,v $


#create the ModelBuilder object
model BasicBuilder -ndm 2 -ndf 2

# build the model 

# add nodes - command: node nodeId xCrd yCrd
node 1   0.0  0.0
node 2 144.0  0.0
node 3 168.0  0.0
node 4  72.0 96.0

# add material - command: uniaxialMaterial <matType> matID <matArgs>
uniaxialMaterial Elastic 1 3000

# add truss elements - command: element truss trussID node1 node2 A matID
#element truss 1 1 4 10.0 1
#element truss 2 2 4 5.0 1
#element truss 3 3 4 5.0 1

# add feap truss elements - command: element fTruss trussID node1 node2 A E
element fTruss 1 1 4 10.0 3000
element fTruss 2 2 4 5.0 3000
element fTruss 3 3 4 5.0 3000

# set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?
fix 1 1 1 
fix 2 1 1
fix 3 1 1

pattern Plain 1 "Linear" {
    # apply the load - command: load nodeID xForce yForce
    load 4 100 -50
}


# build the components for the analysis object
system BandSPD
constraints Plain
integrator LoadControl 1.0 1 1.0 1.0
algorithm Linear
numberer RCM

# create the analysis object 
analysis Static 

# create a Recorder object for the nodal displacements at node 4
recorder Node example.out disp -load -nodes 4 -dof 1 2

# perform the analysis
analyze 1

# print the results at node and at all elements
print node 4
print ele
playback 1

