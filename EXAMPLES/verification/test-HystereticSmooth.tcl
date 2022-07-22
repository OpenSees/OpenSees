wipe
model basic -ndm 2 -ndf 2 



# NODES 
node 1 0.0 0.0 
node 2 1.0 0.0 
mass 2 1.0e1 0.0

set ka 5.0;
set kb 0.5;
set alpha 5;
set beta 1.0

set fbar [expr ($ka-$kb)*0.5/$alpha]


# uniaxialMaterial HystereticSmooth 1 $ka $kb $fbar $beta
uniaxialMaterial HystereticSmooth 1 $ka $kb $fbar $beta



# ELEMENTS 
element truss 1  1 2 1.0  1
 
# BOUNDARY CONDITIONS 
fix 1 1 1 
fix 2 0 1 

recorder Node -file SdofU.out -closeOnWrite -node 2 -dof 1 disp;
recorder Node -file SdofF.out -closeOnWrite -node 1 -dof 1 reaction;
recorder Element -file StressElement.out -closeOnWrite -ele 1 material stress
recorder Element -file StrainElement.out -closeOnWrite -ele 1 material strain
recorder Element -file TangentElement.out -closeOnWrite -ele 1 material tangent

pattern Plain 1 Linear {
load 2 1.0 0.0
}

constraints Plain                            
numberer Plain                               
system BandGeneral                           
test EnergyIncr 1.0e-3 200 ;                
algorithm  KrylovNewton;		 	 # use Kyrlow-Newton algorithm 
#integrator Newmark 0.5 0.25 ;	     # determine the next time step for an analysis
#analysis Transient;					 # define type of analysis: time-dependent                            




integrator DisplacementControl 2 1  0.01   
analysis Static                              
analyze 100

integrator DisplacementControl 2 1  -0.01
analyze 200

integrator DisplacementControl 2 1  0.01
analyze 200