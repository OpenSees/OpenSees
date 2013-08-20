set P 1.0
set L 20.0
set R 1.0
set E 1000.0

set nz  20
set nx  6
set ny  6

set PI [expr 2.0 * asin(1.0)]
set I  [expr $PI*pow((2*$R),4)/64.0] 

puts "PL^3/3EI = [expr $P*pow($L,3)/(3.0*$E*$I)]"

# Create ModelBuilder with 3 dimensions and 6 DOF/node
model Basic -ndm 3 -ndf 3

# create the material
nDMaterial ElasticIsotropic   1   $E   0.25  1.27

set eleArgs "1" 
set element bbarBrick

set nn [expr ($nz)*($nx+1)*($ny+1) + (($nx+1)*($ny+1)+1)/2]
set n1 [expr ($nz)*($nx+1)*($ny+1) +$nx]

# mesh generation
set sqrtR [expr sqrt($R/2.0)]
set cmd "block3D $nx $ny $nz   1 1  $element  $eleArgs {
       1   -$sqrtR     -$sqrtR      0
       2    $sqrtR     -$sqrtR      0
       3    $sqrtR      $sqrtR      0
       4   -$sqrtR      $sqrtR      0 
       5   -$sqrtR     -$sqrtR      $L
       6    $sqrtR     -$sqrtR      $L
       7    $sqrtR      $sqrtR      $L
       8   -$sqrtR      $sqrtR      $L
      13      0         -$R        0
      14      $R         0         0
      15      0          $R        0
      16     -$R         0         0
      18      0         -$R        $L
      19      $R         0         $L
      20      0          $R        $L
      21     -$R         0         $L
      23      0         -$R        [expr $L/2.0]
      24      $R         0         [expr $L/2.0]
      25      0          $R        [expr $L/2.0]
      26     -$R         0         [expr $L/2.0]
    }"

eval $cmd
    
# boundary conditions
fixZ 0.0   1 1 1 

# Constant point load
pattern Plain 1 Linear {
    load $nn  0.0 $P 0.0
}

integrator LoadControl  1.0  
test NormUnbalance     1.0e-10    20     0
algorithm Newton
numberer RCM
constraints Plain 
system ProfileSPD
analysis Static 

analyze 1

set nz [expr 2*$nz]
set ny [expr 2*$ny]
set nx [expr 2*$nx]

print node $nn

recorder display ShakingBeam 100 40 500 500 -wipe
#prp -100 100 120.5
prp 800  0 0
vup 0 0 1
display 1 4 1 

recorder display ShakingBeam2 600 40 500 500 -wipe
#prp -100 100 120.5
prp 800  800 1000
vup 0 0 1
display 1 4 1 


    