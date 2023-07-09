#----------Create Model (2-D 3dof problem)

model basic -ndm 3 -ndf 6

# Create nodes
node 1 0.0 0.0 0.0
node 2 0.0 100.0 0.0


# Fix supports at node
fix 1 1 1 1 1 1 1
fix 2 1 0 1 1 1 1


# Define  material
# uniaxialMaterial BoucWen $matTag $alpha $ko $n $gamma $beta $Ao $deltaA $deltaNu $deltaEta 


#uniaxialMaterial material tag alpha Ko n gama beta A Zetas p Shi deltaShi Lamda tolerance maxNumberIter
#uniaxialMaterial BWBN 1 0.0376 37.76 2 0.5 0.5 1 0.9 1 0.2 0.001 0.05 0.001 1000 

#uniaxialMaterial material tag alpha Ko n gama beta A q Zetas p Shi deltaShi Lamda tolerance maxNumberIter
#uniaxialMaterial BWBN 1 0.0376 37.76 2 0.5 0.5 1 0 0.9 1 0.2 0.001 0.05 0.001 1000 
uniaxialMaterial BWBN 1 0.03 24.877 1 0.1 0.5 1 .4 0.9 1 0.2 0.001 0.05 1e-6 1000 


# Define two node link element

element twoNodeLink 1 1 2 -mat 1 1 1 1 1 1 -dir 1 2 3 4 5 6 -orient 1 0 0 0 1 0


#----------Create Analysis

# Create recorders
recorder Node -file Node2.out -node 2 -dof 2 disp
recorder Element -file Element1.out -ele 1 force


#Cycle 100-3M 1.57 -0.95 1.64 -0.97 1.67 -0.94 3.19 -1.91 3.35 -1.96 3.42 -2.04 7.74 -5.42 7.77 -5.29 7.84 -5.27 -2.92 

pattern Plain 1 Linear {

load 2 0.0 1.0 0.0 0.0 0.0 0.0  

}


constraints Transformation
numberer Plain
integrator LoadControl 1
system BandGeneral
test NormDispIncr 1.0e-6 100
algorithm Newton
analysis Static

integrator DisplacementControl 2 2 0.01
analyze 157
integrator DisplacementControl 2 2 -0.01
analyze 252
integrator DisplacementControl 2 2 0.01
analyze 259
integrator DisplacementControl 2 2 -0.01
analyze 261
integrator DisplacementControl 2 2 0.01
analyze 264
integrator DisplacementControl 2 2 -0.01
analyze 261
integrator DisplacementControl 2 2 0.01
analyze 413
integrator DisplacementControl 2 2 -0.01
analyze 510
integrator DisplacementControl 2 2 0.01
analyze 526
integrator DisplacementControl 2 2 -0.01
analyze 531
integrator DisplacementControl 2 2 0.01
analyze 538
integrator DisplacementControl 2 2 -0.01
analyze 546
integrator DisplacementControl 2 2 0.01
analyze 978
integrator DisplacementControl 2 2 -0.01
analyze 1316
integrator DisplacementControl 2 2 0.01
analyze 1319
integrator DisplacementControl 2 2 -0.01
analyze 1306
integrator DisplacementControl 2 2 0.01
analyze 1313
integrator DisplacementControl 2 2 -0.01
analyze 1311


wipe