# Z. Cheng, UC Davis
# upU test

wipe all

# create the modelbuilder
model BasicBuilder -ndm 3 -ndf 7

# Nodal coordinates
node  1 1.0 1.0 0.0
node  2 0.0 1.0 0.0
node  3 0.0 0.0 0.0
node  4 1.0 0.0 0.0
node  5 1.0 1.0 1.0
node  6 0.0 1.0 1.0
node  7 0.0 0.0 1.0
node  8 1.0 0.0 1.0
node  9 1.0 1.0 2.0
node 10 0.0 1.0 2.0
node 11 0.0 0.0 2.0
node 12 1.0 0.0 2.0
node 13 1.0 1.0 3.0
node 14 0.0 1.0 3.0
node 15 0.0 0.0 3.0
node 16 1.0 0.0 3.0
node 17 1.0 1.0 4.0
node 18 0.0 1.0 4.0
node 19 0.0 0.0 4.0
node 20 1.0 0.0 4.0
node 21 1.0 1.0 5.0
node 22 0.0 1.0 5.0
node 23 0.0 0.0 5.0
node 24 1.0 0.0 5.0
node 25 1.0 1.0 6.0
node 26 0.0 1.0 6.0
node 27 0.0 0.0 6.0
node 28 1.0 0.0 6.0
node 29 1.0 1.0 7.0
node 30 0.0 1.0 7.0
node 31 0.0 0.0 7.0
node 32 1.0 0.0 7.0
node 33 1.0 1.0 8.0
node 34 0.0 1.0 8.0
node 35 0.0 0.0 8.0
node 36 1.0 0.0 8.0
node 37 1.0 1.0 9.0
node 38 0.0 1.0 9.0
node 39 0.0 0.0 9.0
node 40 1.0 0.0 9.0
node 41 1.0 1.0 10.0
node 42 0.0 1.0 10.0
node 43 0.0 0.0 10.0
node 44 1.0 0.0 10.0


# Boundary conditions. 
fix  1  1 1 1 0 1 1 1
fix  2  1 1 1 0 1 1 1
fix  3  1 1 1 0 1 1 1
fix  4  1 1 1 0 1 1 1
fix  5  0 1 0 0 1 1 0
fix  6  0 1 0 0 1 1 0
fix  7  0 1 0 0 1 1 0
fix  8  0 1 0 0 1 1 0
fix  9  0 1 0 0 1 1 0
fix 10  0 1 0 0 1 1 0
fix 11  0 1 0 0 1 1 0
fix 12  0 1 0 0 1 1 0
fix 13  0 1 0 0 1 1 0
fix 14  0 1 0 0 1 1 0
fix 15  0 1 0 0 1 1 0
fix 16  0 1 0 0 1 1 0
fix 17  0 1 0 0 1 1 0
fix 18  0 1 0 0 1 1 0
fix 19  0 1 0 0 1 1 0
fix 20  0 1 0 0 1 1 0
fix 21  0 1 0 0 1 1 0
fix 22  0 1 0 0 1 1 0
fix 23  0 1 0 0 1 1 0
fix 24  0 1 0 0 1 1 0
fix 25  0 1 0 0 1 1 0
fix 26  0 1 0 0 1 1 0
fix 27  0 1 0 0 1 1 0
fix 28  0 1 0 0 1 1 0
fix 29  0 1 0 0 1 1 0
fix 30  0 1 0 0 1 1 0
fix 31  0 1 0 0 1 1 0
fix 32  0 1 0 0 1 1 0
fix 33  0 1 0 0 1 1 0
fix 34  0 1 0 0 1 1 0
fix 35  0 1 0 0 1 1 0
fix 36  0 1 0 0 1 1 0
fix 37  0 1 0 0 1 1 0
fix 38  0 1 0 0 1 1 0
fix 39  0 1 0 0 1 1 0
fix 40  0 1 0 0 1 1 0
fix 41  0 1 0 1 1 1 0
fix 42  0 1 0 1 1 1 0
fix 43  0 1 0 1 1 1 0
fix 44  0 1 0 1 1 1 0

equalDOF  5  6  1 3
equalDOF  5  7  1 3
equalDOF  5  8  1 3
equalDOF  9 10  1 3
equalDOF  9 11  1 3
equalDOF  9 12  1 3
equalDOF 13 14  1 3
equalDOF 13 15  1 3
equalDOF 13 16  1 3
equalDOF 17 18  1 3
equalDOF 17 19  1 3
equalDOF 17 20  1 3
equalDOF 21 22  1 3
equalDOF 21 23  1 3
equalDOF 21 24  1 3
equalDOF 25 26  1 3
equalDOF 25 27  1 3
equalDOF 25 28  1 3
equalDOF 29 30  1 3
equalDOF 29 31  1 3
equalDOF 29 32  1 3
equalDOF 33 34  1 3
equalDOF 33 35  1 3
equalDOF 33 36  1 3
equalDOF 37 38  1 3
equalDOF 37 39  1 3
equalDOF 37 40  1 3
equalDOF 41 42  1 3
equalDOF 41 43  1 3
equalDOF 41 44  1 3




# ##DM04############################################################################################
set e0        0.95
set G0        125.0
set v         0.35
set Pat       100.0
set kc        0.02
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
set rho_s   [expr 2.5/$poro]
set rho_f   1.0
set bulk_s  1.0e12
set bulk_f  2.2e6
set kx      [expr 5e-5/$grvt/$rho_f]
set ky      [expr 5e-5/$grvt/$rho_f]
set kz      [expr 5e-5/$grvt/$rho_f]
set rho     [expr (1.0-$poro)*$rho_s + $poro*$rho_f ]
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

#nDMaterial ElasticIsotropic3D 1 2000 0.3 1.5

#_____________________tag_____8 nodes____matID  bf1___bf2___bf3   poro     alpha     rho_s   rho_f    kx   ky   kz    s_bulk   f_bulk    pressure
element Brick8N_u_p_U  1  5  6  7  8  1  2  3  4   1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0
element Brick8N_u_p_U  2  9 10 11 12  5  6  7  8   1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0
element Brick8N_u_p_U  3  13 14 15 16  9 10 11 12  1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0
element Brick8N_u_p_U  4  17 18 19 20 13 14 15 16  1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0
element Brick8N_u_p_U  5  21 22 23 24 17 18 19 20  1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0
element Brick8N_u_p_U  6  25 26 27 28 21 22 23 24  1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0
element Brick8N_u_p_U  7  29 30 31 32 25 26 27 28  1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0
element Brick8N_u_p_U  8  33 34 35 36 29 30 31 32  1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0
element Brick8N_u_p_U  9  37 38 39 40 33 34 35 36  1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0
element Brick8N_u_p_U 10  41 42 43 44 37 38 39 40  1    0.0   0.0  $Gr  $poro    $alpha    $rho_s  $rho_f   $kx  $ky  $kz   $bulk_s  $bulk_f   0.0

pattern Plain 1 "Constant" {
  eleLoad -ele   1 -type -BrickW
  eleLoad -ele   2 -type -BrickW
  eleLoad -ele   3 -type -BrickW
  eleLoad -ele   4 -type -BrickW  
  eleLoad -ele   5 -type -BrickW
  eleLoad -ele   6 -type -BrickW
  eleLoad -ele   7 -type -BrickW
  eleLoad -ele   8 -type -BrickW  
  eleLoad -ele   9 -type -BrickW
  eleLoad -ele  10 -type -BrickW     
}

set gamma 0.7
integrator Newmark  $gamma  [expr pow($gamma+0.5, 2)/4]
#integrator HHT 0.7 
numberer Plain
constraints Penalty 1.0e12 1.0e12
test NormDispIncr 1.0e-2 30 0
algorithm Newton
system UmfPack
analysis VariableTransient
analyze 300 10  0.1 20

print -node 1 2 3 4 5 6 7 8 9 10 11 12

wipeAnalysis
loadConst -time 0.0


#set Hload [expr 5.0]
#set Sine "Sine 0.2 12.2 0.1"
#pattern Plain 2 $Sine {
# load 1 $Hload 0 0 0 0 0 0
# load 2 $Hload 0 0 0 0 0 0
# load 3 $Hload 0 0 0 0 0 0
# load 4 $Hload 0 0 0 0 0 0
#}

#set Hdisp 0.01
#set Sine "Sine 0 6 0.1 -factor 1.0"
#set Rec "Rectangular 0 2000 -factor 1.0"
#pattern Plain 2 "Linear" {
# sp 25 1 $Hdisp
# sp 26 1 $Hdisp
# sp 27 1 $Hdisp
# sp 28 1 $Hdisp
#}

#===========================================================
recorder Node  -file N01.out  -time -node  25  -dof  1  disp
recorder Node  -file N05.out  -time -node  1  5  9 13 17 21 25 29 33 37 41 -dof  4  disp
recorder Element  -file stress01.out  -time -ele  1 stress
recorder Element  -file stress02.out  -time -ele  2 stress
recorder Element  -file stress05.out  -time -ele  5 stress
recorder Element  -file stress08.out  -time -ele  8 stress

set dt 0.01
#set Gaccel "Series -dt $dt -filePath PSFFTran00.txt -factor [expr $grvt]"
set Gaccel "Series -dt $dt -filePath record01.txt -factor [expr $grvt*20.0]"
pattern UniformExcitation  2  1  -accel $Gaccel

#pattern UniformExcitation  2  1  -accel "Sine 0.1 10.1 0.1 -factor 1.0"

set gamma 0.7
integrator Newmark  $gamma  [expr pow($gamma+0.5, 2)/4]
numberer Plain
constraints Penalty 1.0e16 1.0e16
test NormDispIncr 1.0e-1 30 1
algorithm Newton
system UmfPack
analysis VariableTransient

set dt1 0.01
analyze 3000 $dt1 [expr $dt1/25] $dt1 30

wipe
