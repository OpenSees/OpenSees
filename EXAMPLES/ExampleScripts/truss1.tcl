# written: fmk
# date: 07/99
#
# Simple Truss Example
#
#      +---+---+---+---+
#      |\ /|\ /|\ /|\ /|
#      | X | X | X | X |
#      |/ \|/ \|/ \|/ \|
#      +---+---+---+---+
#      ^       |       0
#     ===      |      ===
#              V  P
#


set DISPLAY OFF

#some variables
set p 5000.0
set E 30000.0
set Fy 60.0
set b 0.1
set A 100.0

set numAnalysisSteps 15

#create the ModelBuilder
model basic -ndm 2 -ndf 2

# define material and section properties
uniaxialMaterial Steel01 3 $Fy $E $b

# define the nodes
node 1    0.0   0.0 
node 2  100.0   0.0 
node 3  200.0   0.0 
node 4  300.0   0.0 
node 5  400.0   0.0 
node 6    0.0 100.0
node 7  100.0 100.0
node 8  200.0 100.0
node 9  300.0 100.0
node 10 400.0 100.0


# define bottom chord elements
element truss 1  1 2 $A 3
element truss 2  2 3 $A 3
element truss 3  3 4 $A 3
element truss 4  4 5 $A 3

# define top chord elements
element truss 5  6 7  $A 3
element truss 6  7 8  $A 3
element truss 7  8 9  $A 3
element truss 8  9 10 $A 3

# define diagonal elements
element truss 9   1 7  $A 3
element truss 10  2 8  $A 3
element truss 11  3 9  $A 3
element truss 12  4 10 $A 3
element truss 13  2 6  $A 3
element truss 14  3 7  $A 3
element truss 15  4 8  $A 3
element truss 16  5 9  $A 3

#define the vertical elements
element truss 17 1  6 $A 3
element truss 18 2  7 $A 3
element truss 19 3  8 $A 3
element truss 20 4  9 $A 3
element truss 21 5 10 $A 3

# fix the left node, right node on rollar
fix 1 1 1 
fix 5 0 1

# create a LoadPattern with a Linear time series
pattern Plain 1 Linear {
    # add a load
    load 3 0.0 -$p
}


#create the recorder
recorder Node -load -node 3 -dof 2 -file Node.out disp
recorder Element -load -eleRange 1 16 -file ElementForce.out axialForce
recorder Element -load -eleRange 1 16 -file ElementStress.out material stress

if {$DISPLAY == "ON"} {
    recorder plot Node.out Node2Disp 50 350 200 200 -columns 2 1
}

# create the SOE, ConstraintHandler, Integrator, Algorithm and Numberer
system SparseGeneral
constraints Plain
integrator LoadControl [expr 1.0/$numAnalysisSteps]
test NormDispIncr 1.0e-8 10 0
algorithm Newton
numberer RCM

# create the Analysis
analysis Static 

if {$DISPLAY == "ON"} {
    # create the display
    recorder display g3 10 10 800 200 -wipe
    prp 20 5.0 100.0
    vup 0 1 0
    viewWindow -30 30 -20 20
    display 1 0 4
}

#analyze the structure
analyze $numAnalysisSteps

puts "Vertical Deflection is [expr -[nodeDisp 3 2]]"





