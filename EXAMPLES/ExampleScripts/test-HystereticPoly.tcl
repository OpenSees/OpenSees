set Time1 [clock seconds]
wipe
model basic -ndm 2 -ndf 2 



# NODES 
node 1 0.0 0.0 
node 2 1.0 0.0 
mass 2 1.0e2 0.0

#source material.tcl
#uniaxialMaterial HystereticPoly $matTag $ka $kb $alpha $beta1 $beta2 <$delta>
uniaxialMaterial HystereticPoly 1  100.0 10.0 20.0 -10.0 10.0 1.0e-8



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

#timeSeries Path 1 -dt 0.01 -filePath EmiliaAccelerogram.th

timeSeries Trig 1 0.0 4.0 0.5 -factor 2.0e0

pattern UniformExcitation 1 1 -accel 1

constraints Plain                            
numberer Plain                               
system BandGeneral                           
test EnergyIncr 1.0e-3 200 ;                
algorithm  KrylovNewton;		 	 # use Kyrlow-Newton algorithm 


integrator DisplacementControl 2 1  0.01   
analysis Static                              

analyze 100

integrator DisplacementControl 2 1  -0.01
analyze 200

integrator DisplacementControl 2 1  0.01
analyze 200