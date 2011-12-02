# Z. Cheng, UC Davis
# upU test

wipe all

# create the modelbuilder
model BasicBuilder -ndm 3 -ndf 7

# Nodal coordinates
node  1 1.0 1.0 1.0
node  2 0.0 1.0 1.0
node  3 0.0 0.0 1.0
node  4 1.0 0.0 1.0
node  5 1.0 1.0 0.0
node  6 0.0 1.0 0.0
node  7 0.0 0.0 0.0
node  8 1.0 0.0 0.0

# Boundary conditions. 
fix  1  0 1 0 1 1 1 0
fix  2  0 1 0 1 1 1 0
fix  3  0 1 0 1 1 1 0
fix  4  0 1 0 1 1 1 0
fix  5  1 1 1 0 1 1 1
fix  6  1 1 1 0 1 1 1
fix  7  1 1 1 0 1 1 1
fix  8  1 1 1 0 1 1 1

equalDOF 1 2  1 3 7
equalDOF 1 3  1 3 7
equalDOF 1 4  1 3 7

##DM04############################################################################################
set e0        0.96
set G0        125.0
set v         0.35
set Pat       100.0
set kc        0.01
set M_cal     1.25
set c         0.8
set lambda_c  0.019
set xi        0.7
set er        0.934
set m_low     0.01
set h0        7.05
set ch        0.968
set nb        1.1
set A0        0.704
set nd        3.5
set z_max     4.0
set cz        600.0
set zZ        "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
set p         0.0
set initS     "$p 0.0 0.0 0.0 $p 0.0 0.0 0.0 $p"

set s1       [expr -(2.0*$v-1.0)/(1.0+$v)]
set s3       [expr -(2.0-4.0*$v)/(1.0+$v)]
set zS       "$s1 0.0 0.0 0.0 $s1 0.0 0.0 0.0 $s3"

## Sepcify some constants
set grvt    9.8
set Gr      [expr -$grvt];
set poro    [expr $e0/(1.0+$e0)]
set alpha   1.0
set rho_s   2.4
set rho_f   1.0
set bulk_s  1.0e12
set bulk_f  2.2e6
set kx      [expr 2e-5/$grvt/$rho_f]
set ky      [expr 2e-5/$grvt/$rho_f]
set kz      [expr 2e-5/$grvt/$rho_f]
#set rho     [expr (1.0-$poro)*$rho_s + $poro*$rho_f ]
set rho     $rho_s
     
##-------1    2   3   4  5    6   7      8  9         10  11  12     13  14  15  16  17  18     19
set mc  "$rho $e0 $G0 $v $Pat $kc $M_cal $c $lambda_c $xi $er $m_low $h0 $ch $nb $A0 $nd $z_max $cz"
set it  "$zS $zZ"
set mp  "MaterialConstant 19 $mc InternalTensor 2 $it"
set el  "Dafalias-Manzari  3 4 5 6 2 $initS"
set yf  "Dafalias-Manzari  0 12 2 1"
set pf  "Dafalias-Manzari  0 2 0 11 0 9 0 10 0 5 0 12 0 7 0 8 0 16 0 17 2 1 2 2"
set te1 "Dafalias-Manzari  2 11 9 10 5 12 7 8 15 13 14 3 1 2"
set te2 "Dafalias-Manzari-fabric  12 19 18 1 2"
set te  "$te1 $te2"

nDMaterial NewTemplate3Dep 1 -MaterialParameter $mp \
                             -ElasticState $el \
                             -YieldFunction $yf \
                             -PlasticFlow $pf \
                             -TensorEvolution $te 



#_____________________tag_____8 nodes____matID  bf1___bf2___bf3   poro     alpha     rho_s   rho_f    kx   ky   kz    s_bulk   f_bulk    pressure
element Brick8N_u_p_U  1  1 2 3 4 5 6 7 8  1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0


pattern Plain 1 "Rectangular 0 1000" {
  eleLoad -ele   1 -type -BrickW
}

recorder Element  -file stress01.out  -time -ele  1 stresses
recorder Node  -file NN05.out  -time -node  5  -dof 4      disp

set gamma 0.6
integrator Newmark  0.6 0.3025
numberer Plain
constraints Penalty 1.0e18 1.0e18
test NormDispIncr 1.0e-3 30 0
algorithm KrylovNewton
system UmfPack
analysis VariableTransient
analyze 200 2

print -node 1 2 3 4 5 6 7 8

wipeAnalysis
loadConst -time 0.0


#===========================================================
recorder Node  -file N01.out  -time -node  1  -dof 1      disp
recorder Node  -file N05.out  -time -node  5  -dof 4      disp
recorder Element  -file stress01.out  -time -ele  1 pq

set dt 0.01

set Gaccel "Series -dt $dt -filePath record01.txt -factor [expr $grvt*10.0]"
pattern UniformExcitation  2  1  -accel $Gaccel


set gamma 0.9
integrator Newmark  $gamma  [expr pow($gamma+0.5, 2)/4]
numberer Plain
constraints Penalty 1.0e20 1.0e20
test NormDispIncr 1.0e-3 30 1
algorithm Newton
system UmfPack
analysis VariableTransient

set dt1 0.005
analyze 10000 $dt1 [expr $dt1/25] $dt1 30

wipe
