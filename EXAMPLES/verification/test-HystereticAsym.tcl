
wipe
model basic -ndm 2 -ndf 2 

# NODES 
node 1 0.0 0.0 
node 2 1.0 0.0 
mass 2 1.0e1 0.0

set ka 5.0;
set kb 0.5;
set alpha 5;
set beta1 1.5
set beta2 1.0
set gamma 0.3

set fo [expr ($ka-$kb)*0.5/$alpha]
puts "fo= $fo"


# uniaxialMaterial HystereticSmooth 1 $ka $kb $fbar $beta
uniaxialMaterial HystereticAsym 1 $ka $kb $fo $beta1 $beta2 $gamma




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

# displacement control analysis
constraints Plain                            
numberer Plain                               
system BandGeneral                           
test EnergyIncr 1.0e-3 200 ;                
algorithm  KrylovNewton;		 	 # use Kyrlow-Newton algorithm                      


integrator DisplacementControl 2 1  0.001   
analysis Static                              
analyze 1000

integrator DisplacementControl 2 1  -0.001
analyze 2000

integrator DisplacementControl 2 1  0.001
analyze 2000


