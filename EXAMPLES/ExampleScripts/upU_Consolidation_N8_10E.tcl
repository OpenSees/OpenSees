# Z. Cheng, UC Davis
# upU element model
# 1D Consolidation 
# Simple elastic example used for verification

wipe all; 

source Units.tcl;

# Parameters
set grvt    [expr 9.8*$m/$sec/$sec ];
set E 			[expr 2.0e4*$kPa ];
set v       0.2;
set poro    0.4;
set alpha   1.0;
set rho_s   [expr 2.0e3*$kg/pow($m,3) ]
set rho_f   [expr 1.0e3*$kg/pow($m,3) ]
set bulk_s  [expr 1.0e20*$Pa ]
set bulk_f  [expr 2.2e9*$Pa ]

# Permeability 
# Might be factored by, for example 'hour' (3600, then real k = k/3600), 
# then the time unit should be considered as 'hour', and not 'second'
set kx      [expr 3.6e-4*$m/$sec/($grvt*$rho_f) ]	
set ky      [expr 3.6e-4*$m/$sec/($grvt*$rho_f) ]	
set kz      [expr 3.6e-4*$m/$sec/($grvt*$rho_f) ]	 

# create the modelbuilder
model BasicBuilder -ndm 3 -ndf 7

# Nodal coordinates
set ht  [expr 1*$m]
set wd  [expr 1*$m]
set Nlayer 2
node   1  [expr $wd]  [expr $wd]  [expr 0*$ht]
node   2         0.0  [expr $wd]  [expr 0*$ht]
node   3         0.0         0.0  [expr 0*$ht]
node   4  [expr $wd]  [expr 0.0]  [expr 0*$ht]
node   5  [expr $wd]  [expr $wd]  [expr 1*$ht]
node   6         0.0  [expr $wd]  [expr 1*$ht]
node   7         0.0         0.0  [expr 1*$ht]
node   8  [expr $wd]  [expr 0.0]  [expr 1*$ht]
node   9  [expr $wd]  [expr $wd]  [expr 2*$ht]
node  10         0.0  [expr $wd]  [expr 2*$ht]
node  11         0.0         0.0  [expr 2*$ht]
node  12  [expr $wd]  [expr 0.0]  [expr 2*$ht]
node  13  [expr $wd]  [expr $wd]  [expr 3*$ht]
node  14         0.0  [expr $wd]  [expr 3*$ht]
node  15         0.0         0.0  [expr 3*$ht]
node  16  [expr $wd]  [expr 0.0]  [expr 3*$ht]
node  17  [expr $wd]  [expr $wd]  [expr 4*$ht]
node  18         0.0  [expr $wd]  [expr 4*$ht]
node  19         0.0         0.0  [expr 4*$ht]
node  20  [expr $wd]  [expr 0.0]  [expr 4*$ht]
node  21  [expr $wd]  [expr $wd]  [expr 5*$ht]
node  22         0.0  [expr $wd]  [expr 5*$ht]
node  23         0.0         0.0  [expr 5*$ht]
node  24  [expr $wd]  [expr 0.0]  [expr 5*$ht]
node  25  [expr $wd]  [expr $wd]  [expr 6*$ht]
node  26         0.0  [expr $wd]  [expr 6*$ht]
node  27         0.0         0.0  [expr 6*$ht]
node  28  [expr $wd]  [expr 0.0]  [expr 6*$ht]
node  29  [expr $wd]  [expr $wd]  [expr 7*$ht]
node  30         0.0  [expr $wd]  [expr 7*$ht]
node  31         0.0         0.0  [expr 7*$ht]
node  32  [expr $wd]  [expr 0.0]  [expr 7*$ht]
node  33  [expr $wd]  [expr $wd]  [expr 8*$ht]
node  34         0.0  [expr $wd]  [expr 8*$ht]
node  35         0.0         0.0  [expr 8*$ht]
node  36  [expr $wd]  [expr 0.0]  [expr 8*$ht]
node  37  [expr $wd]  [expr $wd]  [expr 9*$ht]
node  38         0.0  [expr $wd]  [expr 9*$ht]
node  39         0.0         0.0  [expr 9*$ht]
node  40  [expr $wd]  [expr 0.0]  [expr 9*$ht]
node  41  [expr $wd]  [expr $wd]  [expr 10*$ht]
node  42         0.0  [expr $wd]  [expr 10*$ht]
node  43         0.0         0.0  [expr 10*$ht]
node  44  [expr $wd]  [expr 0.0]  [expr 10*$ht]
    	   
# Boundary Conditions
fix   1  1 1 1 0 1 1 1
fix   2  1 1 1 0 1 1 1
fix   3  1 1 1 0 1 1 1
fix   4  1 1 1 0 1 1 1
fix   5  1 1 0 0 1 1 0
fix   6  1 1 0 0 1 1 0
fix   7  1 1 0 0 1 1 0
fix   8  1 1 0 0 1 1 0
fix   9  1 1 0 0 1 1 0
fix  10  1 1 0 0 1 1 0
fix  11  1 1 0 0 1 1 0
fix  12  1 1 0 0 1 1 0
fix  13  1 1 0 0 1 1 0
fix  14  1 1 0 0 1 1 0
fix  15  1 1 0 0 1 1 0
fix  16  1 1 0 0 1 1 0
fix  17  1 1 0 0 1 1 0
fix  18  1 1 0 0 1 1 0
fix  19  1 1 0 0 1 1 0
fix  20  1 1 0 0 1 1 0
fix  21  1 1 0 0 1 1 0
fix  22  1 1 0 0 1 1 0
fix  23  1 1 0 0 1 1 0
fix  24  1 1 0 0 1 1 0
fix  25  1 1 0 0 1 1 0
fix  26  1 1 0 0 1 1 0
fix  27  1 1 0 0 1 1 0
fix  28  1 1 0 0 1 1 0
fix  29  1 1 0 0 1 1 0
fix  30  1 1 0 0 1 1 0
fix  31  1 1 0 0 1 1 0
fix  32  1 1 0 0 1 1 0
fix  33  1 1 0 0 1 1 0
fix  34  1 1 0 0 1 1 0
fix  35  1 1 0 0 1 1 0
fix  36  1 1 0 0 1 1 0
fix  37  1 1 0 0 1 1 0
fix  38  1 1 0 0 1 1 0
fix  39  1 1 0 0 1 1 0
fix  40  1 1 0 0 1 1 0
fix  41  1 1 0 1 1 1 0
fix  42  1 1 0 1 1 1 0
fix  43  1 1 0 1 1 1 0
fix  44  1 1 0 1 1 1 0

# Material model
nDMaterial ElasticIsotropic3D 1 $E $v $rho_s

# Element model
set Nelement $Nlayer
#_____________________tag______________8 nodes____________matID bf1__bf2__bf3  poro  alpha   rho_s  rho_f   kx  ky  kz   s_bulk   f_bulk   pressure
element Brick8N_u_p_U  1  5  6  7  8  1  2  3  4  1   0.0  0.0  0.0  $poro $alpha  $rho_s $rho_f  $kx $ky $kz  $bulk_s  $bulk_f  0.0
element Brick8N_u_p_U  2  9 10 11 12  5  6  7  8  1   0.0  0.0  0.0  $poro $alpha  $rho_s $rho_f  $kx $ky $kz  $bulk_s  $bulk_f  0.0
element Brick8N_u_p_U  3 13 14 15 16  9 10 11 12  1   0.0  0.0  0.0  $poro $alpha  $rho_s $rho_f  $kx $ky $kz  $bulk_s  $bulk_f  0.0
element Brick8N_u_p_U  4 17 18 19 20 13 14 15 16  1   0.0  0.0  0.0  $poro $alpha  $rho_s $rho_f  $kx $ky $kz  $bulk_s  $bulk_f  0.0
element Brick8N_u_p_U  5 21 22 23 24 17 18 19 20  1   0.0  0.0  0.0  $poro $alpha  $rho_s $rho_f  $kx $ky $kz  $bulk_s  $bulk_f  0.0
element Brick8N_u_p_U  6 25 26 27 28 21 22 23 24  1   0.0  0.0  0.0  $poro $alpha  $rho_s $rho_f  $kx $ky $kz  $bulk_s  $bulk_f  0.0
element Brick8N_u_p_U  7 29 30 31 32 25 26 27 28  1   0.0  0.0  0.0  $poro $alpha  $rho_s $rho_f  $kx $ky $kz  $bulk_s  $bulk_f  0.0
element Brick8N_u_p_U  8 33 34 35 36 29 30 31 32  1   0.0  0.0  0.0  $poro $alpha  $rho_s $rho_f  $kx $ky $kz  $bulk_s  $bulk_f  0.0
element Brick8N_u_p_U  9 37 38 39 40 33 34 35 36  1   0.0  0.0  0.0  $poro $alpha  $rho_s $rho_f  $kx $ky $kz  $bulk_s  $bulk_f  0.0
element Brick8N_u_p_U 10 41 42 43 44 37 38 39 40  1   0.0  0.0  0.0  $poro $alpha  $rho_s $rho_f  $kx $ky $kz  $bulk_s  $bulk_f  0.0
																								 

recorder  Node  -file PorePressure.out  -time -node 1 5 9 13 17 21 25 29 33 37 41 -dof 4 disp
recorder  Node  -file VerticalDisp.out  -time -node 1 5 9 13 17 21 25 29 33 37 41 -dof 3 disp
recorder  Element   -file Stresses.out  -time -ele  1 stresses

# Load Excess Load
set q2    [expr (100.0/4)*$kPa]

set q2tl  [expr -$q2]
set rec   "Rectangular 0.0 1000.0"  
pattern Plain 1 $rec {
load  41  0.0  0.0  $q2tl  0.0  0.0  0.0  0.0 
load  42  0.0  0.0  $q2tl  0.0  0.0  0.0  0.0 
load  43  0.0  0.0  $q2tl  0.0  0.0  0.0  0.0 
load  44  0.0  0.0  $q2tl  0.0  0.0  0.0  0.0 
}

integrator Newmark  0.6  0.3025
numberer RCM
constraints Plain
test NormDispIncr 1.0e-6 10 1
algorithm Newton
system UmfPack
analysis Transient

set     dt    0.1
analyze 3600  $dt
