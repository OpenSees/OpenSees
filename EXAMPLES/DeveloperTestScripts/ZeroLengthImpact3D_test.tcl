
# Prof. Arash E. Zaghi, Majid Cashany, @ UConn
#MAR,02,2012 ADRiK
#units kip, in,


# two blocks on top of each other:

#     y      (z=0.0)                        y     (z=10.0)     
#     |                                     |                 
#     31--------32                          34--------35      
#     |         |                           |         |       
#     21--------22                          24--------25
#
#     11--------12                          14--------15      
#     |         |                           |         |       
#     01--------02---->x                    04--------05---->x

    
#  x      (y=0.0)         x      (y=10.0)         x      (y=10.0)         x      (y=20.0)
#  |                      |                      |                      |
#  02--------05           12--------15           22--------25           32--------35
#  | primary1 |            | primary1 |            |  secondary1 |            |  secondary1 | 
#  01--------04---> z     11--------14---> z     21--------24---> z     31--------34---> z





wipe;
model BasicBuilder -ndm 3 -ndf 3;





# for plane of ( z=0.0) : 
node 01 0.0 0.0 0.0; 
node 02 10.0 0.0 0.0; 
node 11 0.0 10.0 0.0; 
node 12 10.0 10.0 0.0; 
node 21 0.0 10.0 0.0; 
node 22 10.0 10.0 0.0;
node 31 0.0 20.0 0.0; 
node 32 10.0 20.0 0.0;

# for plane of ( z=10.0) : 
node 04 0.0 0.0 10.0; 
node 05 10.0 0.0 10.0; 
node 14 0.0 10.0 10.0; 
node 15 10.0 10.0 10.0;
node 24 0.0 10.0 10.0; 
node 25 10.0 10.0 10.0; 
node 34 0.0 20.0 10.0; 
node 35 10.0 20.0 10.0; 

fix 01 1 1 1;
fix 02 1 1 1;
fix 04 1 1 1;
fix 05 1 1 1;




#nDMaterial ElasticIsotropic $matTag $E $v <$rho>
nDMaterial ElasticIsotropic 1 29000.0 0.3; #steel

#element stdBrick $eleTag $node1 $node2 $node3 $node4 $node5 $node6 $node7 $node8 $matTag
element stdBrick 1 12 11 01 02 15 14 04 05 1;#primary1
element stdBrick 2 32 31 21 22 35 34 24 25 1;#secondary1






################################################## implementation of new element: 
set direction 2; # direction of normal of contact surface
set initGap 0.5; # initial gap
set frictionRatio 0.1; # friction ratio 
set Kt 1.0e5; # penalty stiffness for tangential directions 
set Kn 1.0e5; # penalty stiffness for normal direction 
set Kn2 [expr $Kn * 0.1]; # penalty stiffness after yielding, based on Hertz impact model 
set Delta_y 0.1; # yield displacement based on Hertz impact model  
set cohesion 0.0; # cohesion

element zeroLengthImpact3D 2111 21 11 $direction $initGap $frictionRatio $Kt $Kn $Kn2 $Delta_y $cohesion;
element zeroLengthImpact3D 2212 22 12 $direction $initGap $frictionRatio $Kt $Kn $Kn2 $Delta_y $cohesion;
element zeroLengthImpact3D 2414 24 14 $direction $initGap $frictionRatio $Kt $Kn $Kn2 $Delta_y $cohesion;
element zeroLengthImpact3D 2515 25 15 $direction $initGap $frictionRatio $Kt $Kn $Kn2 $Delta_y $cohesion;
##################################################






#uniaxialMaterial Elastic $matTag $E;
uniaxialMaterial Elastic 6 1.0e-5;
uniaxialMaterial Elastic 7 1.0e-5;

#element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2 ... <-orient $x1 $x2 $x3 $yp1 $yp2 $yp3> <-doRayleigh $rFlag>
element zeroLength 141 11 21 -mat 6 7 6 -dir 1 2 3; # springs with very low stiffness for convergance of Newton-Raphson method 
element zeroLength 142 12 22 -mat 6 7 6 -dir 1 2 3; # springs with very low stiffness for convergance of Newton-Raphson method 
element zeroLength 144 14 24 -mat 6 7 6 -dir 1 2 3; # springs with very low stiffness for convergance of Newton-Raphson method 
element zeroLength 145 15 25 -mat 6 7 6 -dir 1 2 3; # springs with very low stiffness for convergance of Newton-Raphson method 

puts "model is built";





file mkdir output;
recorder Node -file output/disp2Y.out -node 21 22 24 25 -dof 2 disp;
recorder Node -file output/disp1Y.out -node 11 12 14 15 -dof 2 disp;
recorder Element -file output/ele2212g.out -ele 2212 globalForce;


#recorder display $windowTitle $xLoc $yLoc $xPixels $yPixels <-file $fileName>
recorder display "animation" 10 10 800 800 -wipe;
prp 1.0 1.0 1.0;
vup -1.0 1.0 -1.0;#view-up vector (vup)
vpn 1.0 1.0 1.0;#view-plane normal (vpn)
viewWindow -15 15 -15 15;
display 1 3 1;


wipeAnalysis
#system ProfileSPD;
system BandGeneral;
#system BandSPD;
#system SparseGeneral -piv;
#system SparseGeneral;
#constraints Plain;
constraints Penalty 1.0e5 1.0e5;
test NormDispIncr 1.0e-8 1000;#test NormDispIncr $tol $iter <$pFlag>
numberer RCM;
#algorithm Linear;
#algorithm Newton;
#algorithm NewtonLineSearch 0.8;#algorithm NewtonLineSearch $ratio
#algorithm ModifiedNewton;
#algorithm Newton -initial;
algorithm KrylovNewton;
#algorithm BFGS 10;#algorithm BFGS <$count>; $count=number of iterations within a time step until a new tangent is formed
#algorithm Broyden 10;#algorithm Broyden <$count>; $count=number of iterations within a time step until a new tangent is formed

#timeSeries Path $tag -fileTime $fileTime -filePath $filePath <-factor $cFactor>;
timeSeries Path 1 -fileTime example_cyclic_time.txt -filePath example_cyclic_disp.txt;
#pattern Plain $patternTag $tsTag {...}
pattern Plain 1 1 {
            #sp $nodeTag $DOFtag $DOFvalue
            sp 31 2 1.0;
            sp 32 2 1.0;
            sp 34 2 1.0;
            sp 35 2 1.0;
}

#integrator LoadControl $lambda <$numIter $minLambda $maxLambda>;
set maxLambda [expr 0.01];
integrator LoadControl 0.0 1 $maxLambda $maxLambda;
analysis Static;

for {set i 1} {$i<=[expr int(80.0/$maxLambda)]} {incr i 1} {
 algorithm KrylovNewton;
set ok [analyze 1 ];                                     # actually perform analysis; returns ok=0 if analysis was successful


#nodeDisp $nodeTag <$dof>
set a2 [nodeDisp 25 2];
set a1 [nodeDisp 15 2];
#puts " disp. = [expr $a2-$a1]"
 
#puts "$i"
set success 1.0;
if {$ok != 0} {
                        set success 0.0;
                        puts "FAILURE: a miserable one!"
                        break
            }
}


if {$success == 1.0} {
                        puts "SUCCESS: enjoy!"
            }
            
