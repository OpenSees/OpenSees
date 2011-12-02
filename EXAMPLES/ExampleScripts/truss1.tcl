# written: fmk
# date: 07/99
#
# pushover analyusis of truss
#
#      +---+---+---+---+
#      |\ /|\ /|\ /|\ /|
#      | X | X | X | X |
#      |/ \|/ \|/ \|/ \|
#      +---+---+---+---+
#      ^       ^       0
#     ===      |      ===
#              P
#
#
# $Revision: 1.1.1.1 $
# $Date: 2000-09-15 08:23:09 $
# $Source: /usr/local/cvs/OpenSees/EXAMPLES/ExampleScripts/truss1.tcl,v $


#some variables
set p 10
set numAnalysisSteps 15

#create the ModelBuilder
model basic -ndm 2 -ndf 2

# define material and section properties
uniaxialMaterial Elastic 1 300.0
uniaxialMaterial ElasticPP 2 2700.0 0.002
uniaxialMaterial Parallel 3 1 2

# define the nodes
node 1   0.0 0.0 
node 2  10.0 0.0 
node 3  20.0 0.0 
node 4  30.0 0.0 
node 5  40.0 0.0 
node 6   0.0 10.0
node 7  10.0 10.0
node 8  20.0 10.0
node 9  30.0 10.0
node 10 40.0 10.0


# define bottom chord elements
element truss 1  1 2 10.0 3
element truss 2  2 3 10.0 3
element truss 3  3 4 10.0 3
element truss 4  4 5 10.0 3

# define top chord elements
element truss 5  6 7  10.0 3
element truss 6  7 8  10.0 3
element truss 7  8 9  10.0 3
element truss 8  9 10 10.0 3

# define diagonal elements
element truss 9  1 7 10.0 3
element truss 10 2 8 10.0 3
element truss 11 3 9 10.0 3
element truss 12 4 10 10.0 3
element truss 13 2 6  10.0 3
element truss 14 3 7  10.0 3
element truss 15 4 8  10.0 3
element truss 16 5 9  10.0 3

#define the vertical elements
element truss 17 1 6 10.0 3
element truss 18 2 7 10.0 3
element truss 19 3 8 10.0 3
element truss 20 4 9 10.0 3
element truss 21 5 10 10.0 3

# fix the left node, right node on rollar
fix 1 1 1 
fix 5 0 1

# create a LoadPattern with a Linear time series
pattern Plain 1 Linear {
    # add a load
    load 3 0.0 $p
}

#create the recorder
recorder Node Node.out disp -load -node 3  -dof 2
recorder Element 1 -time -file Element.out axialForce
recorder plot Node.out Node2Disp 50 350 200 200 -columns 2 1

# create the SOE, ConstraintHandler, Integrator, Algorithm and Numberer
system SparseGeneral
#constraints Transformation
constraints Plain
integrator LoadControl 1 1 1 1
test NormDispIncr 1.0e-8 10 0
algorithm Newton
numberer RCM

# create the Analysis
analysis Static 

# create the display
recorder display g3 10 10 800 200
prp 20 5.0 100.0
vrp 20 5.0 0
vup 0 1 0
vpn 0 0 1
viewWindow -30 30 -10 10
plane 0 150
port -1 1 -1 1
projection 0
fill 1
display 1 0 5

#analyze the structure
analyze $numAnalysisSteps


