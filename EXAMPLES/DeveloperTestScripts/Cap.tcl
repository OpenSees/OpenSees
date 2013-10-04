
wipe;						                                                           
model basic -ndm 2 -ndf 2;				                                                   
 

node 1 0.0 0.0;					                                                          
node 2 1.0 0.0;
node 3 1.0  1.0
node 4 0.0  1.0;

fix 1 1 1; 
fix 2 1 1;

#          CapPlasticity $tag $ndm $rho     $G           $K          $X     $D        $W   $R   $lambda $theta $beta    $alpha    $T       $tol
nDMaterial CapPlasticity 1    2    2400000  11721092000  14478996000  1.1032e8 4.6412e-10 0.42 4.43 7.9979e6 0.11 6.3816e-8 2.6614e7 -2.0684e6 1.0e-8

set g   9.8
set rhog [expr -2400000*$g]

#set rhog -1.0e12

# 276.8 TON/M3

element quad	1  1  2  3  4	4	"PlaneStrain"	1	0	0.0 	0	$rhog


recorder Node -file node3.out -time -node 3  -dof 1 2 disp
recorder Node -file node4.out -time -node 4  -dof 1 2 disp
recorder Element -ele 1 -time -file stress1.out -time   material 1 stress 
recorder Element -ele 1 -time -file stress2.out -time   material 2 stress 
recorder Element -ele 1 -time -file stress3.out -time   material 3 stress 
recorder Element -ele 1 -time -file stress4.out -time   material 4 stress 

recorder Element  -time -file strain1.out -time -ele 1  material 1 strain
recorder Element -time -file plasticStrain1.out  -ele 1 material 1 plasticStrain 
recorder Element  -time -file tangent1.out -time  -ele 1 material 1 tangent 

recorder Element  -time -file k1.out  -ele 1 material 1 k 
recorder Element  -time -file k2.out -ele 1  material 2 k 
recorder Element  -time -file k3.out  -ele 1 material 3 k 
recorder Element  -time -file k4.out  -ele 1 material 4 k 



#pattern Plain 1 {Series -time {0.0  1000.0 } -values {0.0  10.0 } } {
 
#load 4 1.0e8  -1.0e7 
#load 3 -1.0e7  -1.0e7 
#load 2 -1.0e7  0.0


#}



constraints Plain;     				                                                              
numberer RCM;					                                                              
system BandGeneral;				                                                             
test NormDispIncr 1.0e-8 25 ; 			                                                              
algorithm Newton;					                                                     
integrator LoadControl 0.1;				                                                     
analysis Static	

analyze 10;					                                                            


# -------µØÕð×÷ÓÃ-----------------------------------------------------------------------

loadConst -time 0.0;

wipeAnalysis

#pattern UniformExcitation    1     1    -accel "Series -factor 0.015 -filePath elcentro.txt -dt 0.01"         
pattern UniformExcitation 1 1  -accel  "Series -factor 1.0 -filePath elcentro.txt -dt 0.01"    



constraints Transformation   				                       
numberer Plain					                       
system BandGeneral				
test NormDispIncr 1.0e-8 10 			
algorithm Newton					
integrator Newmark 0.50 0.25 ;			                                                               
analysis Transient

set startT [clock seconds]

analyze 1400 0.01					                                                      


set endT [clock seconds]

puts "Duration: [expr $endT-$startT] seconds."






