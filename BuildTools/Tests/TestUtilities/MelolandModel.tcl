# # Meloland road overcrossing
# 
# ```
#      ^ y
#      |   104'    |    104'   |
#    _____________________________
#    #~o___________o___________o~## __
#    ##\          | |         /####  | 17'
#    z .\_________|o|________/#####------> x
#    ##############################
#
# Nodes:
#
#     4,0          2          1,5
#     ~@-----------o-----------@~
#                  |
#                  o 3
#                 ###
# ```
# 
# - This model assumes linear horizontal and vertical alignment
#   curves
# - Reinforcement is ignored and gross section properties are used
#   in the modeling of the deck.
# 
# The following parameters must be defined by the
# invoking process:
#
# rotation    translation  column   girder
# - k_long    - k_sect     - I_cz   - I_gz
# - k_vert    - k_plan     - I_cx   - I_gy
# - k_tran    - k_elev     - E_c    - E_g
#
#
#         ^ k_plan   sym          ^ k_vert
# k_sect__^___________|___________|_
#    >> __x___________o___________+ -> k_long
#       k_elev       | |        k_tran
#            ________|o|_______
#
# 
#============================================================
# Parameters
#============================================================
set k_long 439679.0
set k_tran 1e+16
set k_vert 2478202.0
set k_sect 1000000.0
set k_plan 2600000.0
set k_elev 20000000000.0
set E_c    4074000.0
set E_g    4074000.0
set A_g    6476.0
set I_cz   636200.0
set I_cy   636200.0
set I_gz   4525000.0
set I_gy   74980000.0
#set bent_mass [lrepeat 3 1301.3]
#set abut_mass [lrepeat 3  500.0]

set bent_mass [lrepeat 3    0.0]
set abut_mass [lrepeat 3  000.0]

set J_col [expr "$I_cz + $I_cy"];

# reduction in accordance with Caltrans 3.4.4-1
set J_col [expr "0.2*$J_col"]; 
set G_col [expr "$E_c   / (2.0*(1+0.2))"];

# Torsional rigidity of the superstructure is not reduced.
set J_gir [expr "$I_gz + $I_gy"]
set G_gir [expr "$E_g   /  (2.0 * (1 + 0.2))"]

#============================================================
# Initializations
#============================================================
# Create ModelBuilder (with 3 dimensions and 6 DOF/node)
model BasicBuilder -ndm 3 -ndf 6
lassign {1 2 3 4 5 6} dx dy dz rx ry rz


#============================================================
# Materials
#============================================================
uniaxialMaterial Elastic  1 $k_long
uniaxialMaterial Elastic  2 $k_vert
uniaxialMaterial Elastic  3 $k_tran
uniaxialMaterial Elastic  4 $k_sect
uniaxialMaterial Elastic  5 $k_plan
uniaxialMaterial Elastic  6 $k_elev
uniaxialMaterial Elastic  7 $k_long
uniaxialMaterial Elastic  8 $k_vert
uniaxialMaterial Elastic  9 $k_tran
uniaxialMaterial Elastic 10 $k_sect
uniaxialMaterial Elastic 11 $k_plan
uniaxialMaterial Elastic 12 $k_elev

#============================================================
# Assemblage
#============================================================
#        x(l)  y(v)   z(t)       x     y     z  s  p  e
node 0      0 240.453   0 -mass  {*}$abut_mass  0  0  0 ; # abut_1
node 1   2496 240.453   0 -mass  {*}$abut_mass  0  0  0 ; # abut_3
node 2   1248 240.453   0 -mass  {*}$bent_mass  0  0  0 ; # bent_top
node 3   1248      0    0 -mass  0     0     0  0  0  0 ; # bent_bot
node 4      0 240.453   0 -mass  0     0     0  0  0  0 ; # abut_1-copy
node 5   2496 240.453   0 -mass  0     0     0  0  0  0 ; # abut_3-copy
geomTransf Linear 1 0 0 1 
geomTransf Linear 2 0 0 1 -jntOffset 0.0 0.0 0.0 0.0 -36.45278977953934 0.0


#                         t i j      A      E      G      J     Iy      Iz   <T>
element elasticBeamColumn 0 0 2    $A_g   $E_g  $G_gir $J_gir  $I_gy  $I_gz   1  -mass 1.4074228 -cMass
element elasticBeamColumn 1 3 2 2.827e+03 $E_c  $G_col $J_col  $I_cy  $I_cz   2  -mass 0.6145110 -cMass
element elasticBeamColumn 2 2 1    $A_g   $E_g  $G_gir $J_gir  $I_gy  $I_gz   1  -mass 1.4074228 -cMass

element zeroLength 3   4   0 -mat 1 2 3  4  5  6 -dir $dx $dy $dz $rx $ry $rz 
element zeroLength 4   5   1 -mat 7 8 9 10 11 12 -dir $dx $dy $dz $rx $ry $rz 

#============================================================
# Constraints
#============================================================
fix 3 1 1 1 1 1 1
fix 4 1 1 1 1 1 1
fix 5 1 1 1 1 1 1

#============================================================
# Eigenvalue analysis
# This file conducts an eigenvalue analysis using OpenSees 
# and writes the results to a file named `results.out`.
#
# Claudio Perez
# Summer 2021
# OpenSees version 3.3.0
#============================================================
# Constant parameters.
set verbose 0
set PI 3.1415159
set numModes 3
set DOFs     {1 2 3 4 5 6}
set nodeRange {0 2}
# Initialize variables `omega`, `f` and `T` to
# empty lists.
foreach {omega f T} {{} {} {}} {}

for {set k 1} {$k <= $numModes} {incr k} {
  recorder Node -nodeRange {*}$nodeRange -dof {*}$DOFs "eigen $k";
}

set eigenvals [eigen $numModes];

foreach eig $eigenvals {
  lappend omega [expr sqrt($eig)];
  lappend f     [expr sqrt($eig)/(2.0*$PI)];
  lappend T     [expr 1000.0*(2.0*$PI)/sqrt($eig)];
}

# print info to `stdout`.
if $verbose {
  puts "Angular frequency (rad/s): \t$omega\n";
  puts "Frequency (Hz): \t$f\n";
  puts "Periods (sec): \t$T\n";
}

proc is_tran {mode_num} {
  global dx dz verbose
  set bent_top_tran [nodeEigenvector 2 $mode_num $dz]
  set bent_top_long [nodeEigenvector 2 $mode_num $dx]
  if $verbose {puts "eigvect: $bent_top_tran $bent_top_long"}
  return [expr {abs($bent_top_tran) > abs($bent_top_long)}]
}
proc is_plan {mode_num} {
  global dx ry verbose
  set bent_top_plan [nodeEigenvector 2 $mode_num $ry]
  set bent_top_long [nodeEigenvector 2 $mode_num $dx]
  return [expr {abs($bent_top_plan) > abs($bent_top_long)}]
}

# If mode number 1 is transverse, long_index is `1` (2nd mode)
# set tran_index [expr ![is_tran 1]]
set tran_index 0
set long_index [expr  [is_plan 2]+1]

# Print second eigenvector to stdout for later processing;
# when using quoFEM these values can be obtained from the
# `ops.out` files that are left in the working directory.
set mode_1 "[nodeEigenvector 0 1] [nodeEigenvector 2 1] [nodeEigenvector 1 1]"
set mode_2 "[nodeEigenvector 0 2] [nodeEigenvector 2 2] [nodeEigenvector 1 2]"
set mode_3 "[nodeEigenvector 0 3] [nodeEigenvector 2 3] [nodeEigenvector 1 3]"
#puts "[lindex $T 0] [lindex $T 1] [lindex $T 2] $mode_1 $mode_2 $mode_3"
puts "[lindex $T $tran_index] [lindex $T $long_index] 0.0 $mode_1 $mode_2 $mode_3"

# Open output file and write first eigenvalue.
set OutputFile [open results.out w]
puts $OutputFile "[lindex $T $tran_index] [lindex $T $long_index]";
#puts $OutputFile "[lindex $T 0] [lindex $T 1]";
close $OutputFile


