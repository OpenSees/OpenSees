# unit kPa kN m s
wipe

model BasicBuilder -ndm 3 -ndf 3

#nDmaterial Simplified3DJ2  $tag      $G          $K         $sig0       $H_kin     $H_iso 
nDMaterial Simplified3DJ2    2       0.8e8       1.75e8     335000       2.1e6      0 ;     
# nDMaterial Simplified3DJ2    2     1.2e8       1.7e9      248200       1e7       1e7

#         create meshes
node        1      0.00000     0.0000    0.0000
node        2      1.00000     0.0000    0.0000
node        3      1.00000     1.0000    0.0000
node        4      0.00000     1.0000    0.0000

node        5      0.00000     0.0000    1.0000
node        6      1.00000     0.0000    1.0000
node        7      1.00000     1.0000    1.0000
node        8      0.00000     1.0000    1.0000



element bbarBrick      1      1    2   3   4    5   6   7   8   2               
            

fix      1      1      1      1 
fix      2      1      1      1  
fix      3      1      1      1   
fix      4      1      1      1   

recorder Element      -file   strain_5.out        -time          -ele   1    material  5 strain ;			   
recorder Element      -file   stress_5.out        -time          -ele   1    material  5   stress ;
recorder Element      -file   plasticStrainDev_5.out        -time          -ele   1    material  5   plasticStrainDev ;



set N  1e3;
 
pattern Plain 1 "Linear" {


    load 5         [expr  $N/4]        0        0    
    load 6         [expr  $N/4]        0        0    
    load 7         [expr  $N/4]        0        0    
    load 8         [expr  $N/4]        0        0    
}


system BandGeneral
test NormDispIncr 1.0e-8 10 
constraints Transformation
integrator LoadControl 1 1 1 1
algorithm Newton 
numberer RCM
analysis Static
analyze 200
                
 

