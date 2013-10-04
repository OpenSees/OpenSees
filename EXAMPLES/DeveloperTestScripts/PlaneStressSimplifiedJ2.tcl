# unit kPa kN m s
wipe

model BasicBuilder -ndm 2 -ndf 2

#nDmaterial                $tag      $G          $K         $sig0       $H_kin     $H_iso 
nDMaterial PlaneStressSimplifiedJ2    2       0.8e8       1.75e8     335000       2.1e6      0 ;     
 


#         create meshes
node        1      0.00000     0.0000    
node        2      1.00000     0.0000   
node        3      1.00000     1.0000   
node        4      0.00000     1.0000   

 

#element FourNodeQuad eleTag? iNode? jNode? kNode? lNode? thk?   type? matTag
element      quad      1      1        2     3       4      0.1  "PlaneStress"       2               
            

fix      1      1      1       
fix      2      1      1       


recorder Element      -file   strain_1.out        -time          -ele   1    material 2  strain ;			   
recorder Element      -file   stress_1.out        -time          -ele   1    material 2  stress ;
recorder Element      -file   strain33.out         -time          -ele   1    material 2  strain33 ;



set N  1e3;
 
pattern Plain 1 "Linear" {


    load 3         [expr  $N/4]        0            
    load 4         [expr  $N/4]        0             
   
}


system BandGeneral
test NormDispIncr 1.0e-8 10 
constraints Transformation
integrator LoadControl 1 1 1 1
algorithm Newton 
numberer RCM
analysis Static
analyze 500
                
 

