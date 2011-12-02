################################################################################
#                                                                              #
# This program conducts nonlinear dynamic analysis for bridge and soil system  #
#                                                                              #
# By: Yuyi Zhang (UCLA) and Zhaohui Yang (UCSD)                                #
# Date: 06/25/01                                                               #
#                                                                              #
################################################################################


# Units: KN, m, sec


proc addNode {nodeNum xDim yDim nodeOutput} {
    global X Y
    incr nodeNum 1
    node $nodeNum $xDim $yDim 
    set X($nodeNum) $xDim; set Y($nodeNum) $yDim;
    puts $nodeOutput "$nodeNum $xDim $yDim"
    expr $nodeNum
}

wipe

#----------------------------UCLA---------------------------------
set in2m 2.54e-2
set kips2Ton 0.4535929
set g 9.81
set kips2KN [expr $kips2Ton*$g]


set wdc 8.6806e-5;     # Weight density of concrete in the columns, kips/in3      
set wdc [expr $wdc*$kips2Ton/($in2m*$in2m*$in2m)]; # Units: ton
set pierWidth 84.0 ;   # inch
set pierWidth [expr $pierWidth*$in2m]
set pierDepth 48.0 ;   # inch
set pierDepth [expr $pierDepth*$in2m]
set triDepth 21.0  ;   # inch
set triDepth [expr $triDepth*$in2m]
set cover 3.0      ;   # inch
set cover [expr $cover*$in2m]
set As 1.48;           # Area of no. 11 bar in the columns, in2
set As [expr $As*$in2m*$in2m]
set fy 60.0;           # Yield strength of reinforcing steel, ksi
set fy [expr $fy*$kips2KN/($in2m*$in2m)]
set E 29000;           # Young's modulus of reinforcing steel, ksi
set E [expr $E*$kips2KN/($in2m*$in2m)]
set np 5;              # Number of Gauss-Lobato points per beam-column element
set fcCore -5.0;       # f'c of core concrete, ksi
set fcCore [expr $fcCore*$kips2KN/($in2m*$in2m)]
set fcuCore -3.0;      # f'cu of core concrete, ksi 
set fcuCore [expr $fcuCore*$kips2KN/($in2m*$in2m)]
set fcCover -4.0;      # f'c of cover concrete, ksi
set fcCover [expr $fcCover*$kips2KN/($in2m*$in2m)]
set fcuCover 0.0;      # f'cu of cover concrete, ksi
set fcuCover [expr $fcuCover*$kips2KN/($in2m*$in2m)]
set beamSectionArea 1338.0; #in2
set beamSectionArea [expr $beamSectionArea*$in2m*$in2m]
set columnSectionArea [expr $pierWidth*$pierDepth+$triDepth*$pierDepth]
set beamInertia 1.01e6;    #inch4
set beamInertia [expr $beamInertia*$in2m*$in2m*$in2m*$in2m]
set beamE 4030.0 ;     # ksi
set beamE [expr $beamE*$kips2KN/($in2m*$in2m)]


#set colLength1 418.6
#set colLength1 [expr $colLength1*$in2m]
#set colLength2 348.5
#set colLength2 [expr $colLength2*$in2m]
#set colLength3 394.4
#set colLength3 [expr $colLength3*$in2m]
#set colLength4 415.1
#set colLength4 [expr $colLength4*$in2m]
#set colLength5 419.5
#set colLength5 [expr $colLength5*$in2m]
#set colLength6 409.8
#set colLength6 [expr $colLength6*$in2m]
#set colLength7 547.4
#set colLength7 [expr $colLength7*$in2m]
#set colLength8 439.2
#set colLength8 [expr $colLength8*$in2m]

set colLength1 10.2
set colLength2 11.7
set colLength3 12.9
set colLength4 13.5
set colLength5 13.6
set colLength6 13.3
set colLength7 12.5
set colLength8 11.4

set tol 1.0e-3

set nodeNumMaxUCSD 1245
set eleNumMaxUCSD 1122
set matNumMaxUCSD 18
set transfNumMaxUCSD 2

set abutmentLeft 1155
set abutmentRight 1158
set col1BaseNode 1188
set col2BaseNode 1196
set col3BaseNode 1202
set col4BaseNode 1210
set col5BaseNode 1216
set col6BaseNode 1224
set col7BaseNode 1232
set col8BaseNode 1238 
#---------------------------UCLA----------------------------------



#some user defined variables
# 
set massDen  2      ;# soil mass density
set accNam  rrsmv1_1.txt ;# acc. file name
set accMul   9.81    ;# acc. time history multiplier
set accDt   0.005  ;# time step of input acc
set deltaT   0.005   ;# time step for analysis
set numSteps 2     ;# number of time steps
set alpha    0.6    ;# Newmark integration parameter
set crack    1      ;# use cracked concrete section
set z 3.281  ;#conversion from feet to meter
#set thick1 [expr 46./$z]  ;# thickness of soil mesh
set thick1 [expr 20./$z]  ;# thickness of soil mesh
set thick2 [expr 34./$z]  ;# thickness of abutment mesh

set nodeOutput [open node w]  ;# node number and coordinates output
set eleOutput [open ele w]    ;# element number and nodes output

#############################################################
# Define model builder for quads
model basic -ndm 2 -ndf 2

# define material and properties
# 1/5 SP to SP/SM (dense) assume average N60 = 50 => Dr=80%, phi=40 degrees.
# sigma_m = 70 kPa
    nDMaterial PressureDependMultiYield 1 2 1.e5 3.e5 40. .1 70. 1. 0.5 25 \
                                        30.0 0.05 -0.00 0.05 20. 0.04 0 0.0 0.0 0. 101

# 2/6 SP/SM (medium) assume average N60 = 20 => Dr=40%, phi=32 degrees.
    nDMaterial PressureDependMultiYield 2 2 4.e4 1.2e5 32. .1 80. 1. 0.5 25 \
                                        27.0 0.05 -0.05 0.05 10. 0.04 20 0.1 3.0 1. 101

# 3/7 OL to OL/SM assume average N60 = 6 => phi=25 degrees.
    nDMaterial PressureDependMultiYield 3 2 3.e4 1.e5 30. .1 80. 1. 0.5 25 \
                                        30.0 0.02 -0.00 0.00 0. 0.04 0 0.1 3.0 0. 101

# 4/8 OL assume average N60 = 10 => phi=25 degrees.
    nDMaterial PressureDependMultiYield 4 2 3.e4 1.e5 25. .1 80 1 0.5 25 \
                                        25.0 0.02 -0.00 0.05 0. 0.04 0 0.1 3.0 0. 101
    nDMaterial FluidSolidPorous 5 2 1 2.2e6 101  ;#undrained, pressure-depend
    nDMaterial FluidSolidPorous 6 2 2 2.2e6 101  ;#undrained, pressure-depend
    nDMaterial FluidSolidPorous 7 2 3 2.2e6 101  ;#undrained, pressure-depend
    nDMaterial FluidSolidPorous 8 2 4 2.2e6 101  ;#undrained, pressure-depend

# 9/10 CL assume average N60 = 20 => phi=30 degrees, Su=20 kPa
    nDMaterial PressureIndependMultiYield 9 2 2.e4 1.e5 30. .1 80 20. 0.5 25 
    nDMaterial FluidSolidPorous 10 2 9 2.2e6 101  ;#undrained, pressure-independ

# Abutment
    nDMaterial PressureIndependMultiYield 11 2 3.e4 1.e5 30. .1 80 30 0.5 25 
    nDMaterial ElasticIsotropic 12 2000 0.25

    for {set i 1} {$i <=11 } {incr i 1} {
       updateMaterialStage -material $i -stage 0
    }

    array set eleMat {
       1 PlaneStrain
       2 PlaneStrain
       3 PlaneStrain
       4 PlaneStrain
       5 PlaneStrain
       6 PlaneStrain
       7 PlaneStrain
       8 PlaneStrain
       9 PlaneStrain
      10 PlaneStrain
      11 PlaneStrain
      12 PlaneStrain
    }

    set grav1 [expr -9.81*$massDen]  ;#gravity
    set grav2 [expr -9.81*($massDen-1.0)]  ;# buoyant unit weight
    for {set i 1} {$i <=12 } {incr i 1} {
       if {[string compare $eleMat($i) FluidSolidPorous]==0} {
          set gravY($i) $grav2;
       } else {
          set gravY($i) $grav1;
       }
    }

# 5 6 7 10 8 11
array set matter {
  1  5 
  2  6
  3  7
  4 10  
  5  8
  6 11
}

# sets the number of blocks (super elements) in two direction.
set numXBlk 16
set numYBlk 8

# sets the list of number of elements in two directions
set numXEleList [list 6 12 4 4 2 5 9 6 6 3 3 3 3 6 12 6]
if {[llength $numXEleList] != $numXBlk} {
   error "Insufficient numbers in numXEleList"
} 
for {set i 0} {$i < $numXBlk } {incr i 1} {
   set numXEle($i) [lindex $numXEleList $i]
}

set numYEleList [list 3 2 1 1 1 1 1 1]
if {[llength $numYEleList] != $numYBlk} {
   error "Insufficient numbers in numYEleList"
} 
for {set i 0} {$i < $numYBlk } {incr i 1} {
   set numYEle($i) [lindex $numYEleList $i]
}
unset numXEleList numYEleList

# coordinates of control nodes 
set x0 0.; set x1 [expr $x0+36.]; set x2 [expr $x1+72.]; set x3 [expr $x2+24.5]; 
set x4 [expr $x3+24.5]; set x5 [expr $x4+12.2]; set x6 [expr $x5+30.5]; 
set x7 [expr $x6+54.9]; set x8 [expr $x7+36.6]; set x9 [expr $x8+36.6]; 
set x10 [expr $x9+18.3]; set x11 [expr $x10+18.3]; set x12 [expr $x11+18.3]; 
set x13 [expr $x12+18.3]; set x14 [expr $x13+36.]; set x15 [expr $x14+72.]; 
set x16 [expr $x15+36.];

set xLeftMost $x0; set xRightMost $x16;

set XList(0) [list $x0 $x0 $x0 $x0 $x0 $x0 $x0 $x0 $x0]
set XList(1) [list $x1 $x1 $x1 $x1 $x1 $x1 $x1 $x1 $x1]
set XList(2) [list $x2 $x2 $x2 $x2 $x2 $x2 $x2 $x2 $x2]
set XList(3) [list $x3 $x3 $x3 $x3 $x3 $x3 $x3 $x3 $x3]
set XList(4) [list $x4 $x4 $x4 $x4 $x4 $x4 $x4 $x4 $x4]
set XList(5) [list $x5 $x5 $x5 $x5 $x5 $x5 $x5 $x5 $x5]
set XList(6) [list $x6 $x6 $x6 $x6 $x6 $x6 $x6 $x6 $x6]
set XList(7) [list $x7 $x7 $x7 $x7 $x7 $x7 $x7 $x7 $x7]
set XList(8) [list $x8 $x8 $x8 $x8 $x8 $x8 $x8 $x8 $x8]
set XList(9) [list $x9 $x9 $x9 $x9 $x9 $x9 $x9 $x9 $x9]
set XList(10) [list $x10 $x10 $x10 $x10 $x10 $x10 $x10 $x10 $x10]
set XList(11) [list $x11 $x11 $x11 $x11 $x11 $x11 $x11 $x11 $x11]
set XList(12) [list $x12 $x12 $x12 $x12 $x12 $x12 $x12 $x12 $x12]
set XList(13) [list $x13 $x13 $x13 $x13 $x13 $x13 $x13 $x13 $x13]
set XList(14) [list $x14 $x14 $x14 $x14 $x14 $x14 $x14 $x14 $x14]
set XList(15) [list $x15 $x15 $x15 $x15 $x15 $x15 $x15 $x15 $x15]
set XList(16) [list $x16 $x16 $x16 $x16 $x16 $x16 $x16 $x16 $x16]

set yb [expr 120./$z]  ;#base elevation -120 ft
set YList(0) [list 0. [expr $yb-80./$z] [expr $yb-63./$z] [expr $yb-55./$z] [expr $yb-45./$z] [expr $yb-32./$z] [expr $yb-19./$z] [expr $yb-8./$z] $yb]
set YList(1) [list 0. [expr $yb-80./$z] [expr $yb-63./$z] [expr $yb-55./$z] [expr $yb-45./$z] [expr $yb-32./$z] [expr $yb-19./$z] [expr $yb-8./$z] $yb]
set YList(2) [list 0. [expr $yb-80./$z] [expr $yb-58./$z] [expr $yb-47./$z] [expr $yb-35./$z] [expr $yb-27./$z] [expr $yb-18./$z] [expr $yb-8./$z] $yb]
set YList(3) [list 0. [expr $yb-80./$z] [expr $yb-60./$z] [expr $yb-50./$z] [expr $yb-40./$z] [expr $yb-32./$z] [expr $yb-22./$z] [expr $yb-10./$z] [expr $yb-5./$z]]
set YList(4) [list 0. [expr $yb-80./$z] [expr $yb-63./$z] [expr $yb-55./$z] [expr $yb-45./$z] [expr $yb-35./$z] [expr $yb-25./$z] [expr $yb-12./$z] [expr $yb-7./$z]]
set YList(5) [list 0. [expr $yb-80./$z] [expr $yb-57./$z] [expr $yb-45./$z] [expr $yb-40./$z] [expr $yb-30./$z] [expr $yb-20./$z] [expr $yb-15./$z] [expr $yb-10./$z]]
set YList(6) [list 0. [expr $yb-80./$z] [expr $yb-59./$z] [expr $yb-48./$z] [expr $yb-45./$z] [expr $yb-38./$z] [expr $yb-30./$z] [expr $yb-23./$z] [expr $yb-17./$z]]
set YList(7) [list 0. [expr $yb-80./$z] [expr $yb-59./$z] [expr $yb-48./$z] [expr $yb-45./$z] [expr $yb-38./$z] [expr $yb-32./$z] [expr $yb-25./$z] [expr $yb-23./$z]]
set YList(8) [list 0. [expr $yb-80./$z] [expr $yb-59./$z] [expr $yb-48./$z] [expr $yb-45./$z] [expr $yb-38./$z] [expr $yb-30./$z] [expr $yb-23./$z] [expr $yb-21./$z]]
set YList(9) [list 0. [expr $yb-80./$z] [expr $yb-60./$z] [expr $yb-50./$z] [expr $yb-45./$z] [expr $yb-38./$z] [expr $yb-31./$z] [expr $yb-24./$z] [expr $yb-18./$z]]
set YList(10) [list 0. [expr $yb-82./$z] [expr $yb-60./$z] [expr $yb-50./$z] [expr $yb-40./$z] [expr $yb-33./$z] [expr $yb-26./$z] [expr $yb-20./$z] [expr $yb-13./$z]]
set YList(11) [list 0. [expr $yb-85./$z] [expr $yb-62./$z] [expr $yb-50./$z] [expr $yb-40./$z] [expr $yb-30./$z] [expr $yb-23./$z] [expr $yb-15./$z] [expr $yb-7./$z]]
set YList(12) [list 0. [expr $yb-86./$z] [expr $yb-62./$z] [expr $yb-50./$z] [expr $yb-40./$z] [expr $yb-30./$z] [expr $yb-20./$z] [expr $yb-10./$z] $yb]
set YList(13) [list 0. [expr $yb-88./$z] [expr $yb-63./$z] [expr $yb-50./$z] [expr $yb-40./$z] [expr $yb-30./$z] [expr $yb-20./$z] [expr $yb-10./$z] $yb]
set YList(14) [list 0. [expr $yb-88./$z] [expr $yb-63./$z] [expr $yb-50./$z] [expr $yb-40./$z] [expr $yb-30./$z] [expr $yb-20./$z] [expr $yb-10./$z] $yb]
set YList(15) [list 0. [expr $yb-88./$z] [expr $yb-63./$z] [expr $yb-50./$z] [expr $yb-40./$z] [expr $yb-30./$z] [expr $yb-20./$z] [expr $yb-10./$z] $yb]
set YList(16) [list 0. [expr $yb-88./$z] [expr $yb-63./$z] [expr $yb-50./$z] [expr $yb-40./$z] [expr $yb-30./$z] [expr $yb-20./$z] [expr $yb-10./$z] $yb]

# sets index for blocks with diagonals 
# 0 => no diagonal; 1 => "slash" diagonal; 2=> "back slash" diagonal
for {set i 0} {$i < $numXBlk} {incr i 1} {
   for {set j 0} {$j < $numYBlk} {incr j 1} {
      set diag($i,$j,0) 0  ;# no diagonal
      set diag($i,$j,1) $matter(1)  ;# material number
      if {$i>=11} { set diag($i,$j,1) $matter(3) }
      if {$j==[expr $numYBlk-1]} { set diag($i,$j,1) $matter(3) }
      if {$j==0} { set diag($i,$j,1) $matter(5) }
      if {$j==1} { set diag($i,$j,1) $matter(1) }
      if { ($j==2) && ($i<14) } { set diag($i,$j,1) $matter(1) }
      if { ($j==2) && ($i>=14) } { set diag($i,$j,1) $matter(3) }
      if { ($j==3) && ($i<=7) } { set diag($i,$j,1) $matter(4) }
      if { ($j==6) && ($i<=4) } { set diag($i,$j,1) $matter(2) }
      if { ($j>=4) && ($j<=5) && ($i<=2) } { set diag($i,$j,1) $matter(2) }
   }
}
  
# diag(?,?,1) material number for lower triangle
# diag(?,?,2) material number for upper triangle
set diag(2,4,0) 1; set diag(2,4,1) $matter(1); set diag(2,4,2) $matter(2);  
set diag(3,5,0) 1; set diag(3,5,1) $matter(1); set diag(3,5,2) $matter(2); 
set diag(4,6,1) $matter(2); 

set diag(9,6,0) 2; set diag(9,6,1) $matter(1);  set diag(9,6,2) $matter(3);  
set diag(10,5,0) 2; set diag(10,5,1) $matter(1);  set diag(10,5,2) $matter(3); 
set diag(10,6,1) 3; 
set diag(11,3,1) 1;
set diag(11,4,0) 2; set diag(11,4,1) $matter(1);  set diag(11,4,2) $matter(3);  
set diag(12,3,0) 2; set diag(12,3,1) $matter(1);  set diag(12,3,2) $matter(3); 
set diag(13,2,0) 2; set diag(13,2,1) $matter(1);  set diag(13,2,2) $matter(3); 

for {set i 0} {$i < $numXBlk} {incr i 1} {
   for {set j 0} {$j < $numYBlk} {incr j 1} {
      if { ($diag($i,$j,0) != 0) && (fmod($numXEle($i),$numYEle($j)) != 0) } {
         error "Block ($i,$j) number of elements in x must be multiple of that in y"
      }
   }
}

# generate nodes and elements
set nodeNum 0
set eleNum 0
set lnode {}
set rnode {}
set mat {};

for {set i 0} {$i < $numYBlk} {incr i 1} {
   set xStart [lindex $XList(0) $i]
   set xEnd [lindex $XList(0) [expr $i+1]]
   set yStart [lindex $YList(0) $i]
   set yEnd [lindex $YList(0) [expr $i+1]]
   set xIncr [expr ($xEnd - $xStart)/$numYEle($i) ]
   set yIncr [expr ($yEnd - $yStart)/$numYEle($i) ]
   for {set j 0} {$j < $numYEle($i)} {incr j 1} {
      set xDim [expr $xStart + $j*$xIncr]
      set yDim [expr $yStart + $j*$yIncr]
      set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
      lappend lnode $nodeNum

      # add one more node at diagonal 
      if { ($diag(0,$i,0)==1) && ($j==0) } {
         set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
         equalDOF [expr $nodeNum-1] $nodeNum 1 2
         lappend lnode $nodeNum
      } 
      if { ($i > 0) && ($diag(0,[expr $i-1],0)==2) && ($j==0) } {
         set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
         equalDOF [expr $nodeNum-1] $nodeNum 1 2
         lappend lnode $nodeNum
      } 
   }
   # last block
   if { $i == [expr $numYBlk-1] } {
      set xDim $xEnd
      set yDim $yEnd
      set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
      lappend lnode $nodeNum
      if { $diag(0,$i,0)==2 } {
         set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
         equalDOF [expr $nodeNum-1] $nodeNum 1 2
         lappend lnode $nodeNum 
      }
   }
}

for {set i 0} {$i < $numXBlk} {incr i 1} {
   for {set k 1} {$k <= $numXEle($i)} {incr k 1} {

      # rearrange lnode
      if { ($i>0) && ($k==1) } {
   
         # remove diagonal-induced equalDOF nodes
         set nCount 0 
         for {set j 0} {$j < $numYBlk} {incr j 1} {
            set im1 [expr $i-1]
            set newCount [expr $nCount+$numYEle($j)]
            if { $diag($im1,$j,0)==2 } {
               set lnode [lreplace $lnode $nCount $nCount]
            } elseif { $diag($im1,$j,0)==1 } {
               set lnode [lreplace $lnode $newCount $newCount]
            }
            set nCount $newCount
         } 
      
         # insert equalDOF nodes due to diagonal
         set nCount 0 
         for {set j 0} {$j < $numYBlk} {incr j 1} {
            set newCount [expr $nCount+$numYEle($j)]
            if { $diag($i,$j,0)==1 } {
               set nod [lindex $lnode $nCount]
               set nodeNum [addNode $nodeNum $X($nod) $Y($nod) $nodeOutput] 
               equalDOF $nod $nodeNum 1 2
               set lnode [linsert $lnode $nCount $nodeNum]
               incr newCount 1
            } elseif { $diag($i,$j,0)==2 } {
               set nod [lindex $lnode $newCount]
               set nodeNum [addNode $nodeNum $X($nod) $Y($nod) $nodeOutput] 
               equalDOF $nod $nodeNum 1 2
               set lnode [linsert $lnode $newCount $nodeNum]
               incr newCount 1
            }
            set nCount $newCount
         } 
      }

      for {set j 0} {$j < $numYBlk} {incr j 1} {
         set jp1 [expr $j+1]
         set ip1 [expr $i+1]
         set xLeftStart [lindex $XList($i) $j]
         set xLeftEnd [lindex $XList($i) $jp1]
         set yLeftStart [lindex $YList($i) $j]
         set yLeftEnd [lindex $YList($i) $jp1]
         set xRightStart [lindex $XList($ip1) $j]
         set xRightEnd [lindex $XList($ip1) $jp1]
         set yRightStart [lindex $YList($ip1) $j]
         set yRightEnd [lindex $YList($ip1) $jp1]
         set xStart [expr $xLeftStart + $k*($xRightStart-$xLeftStart)/$numXEle($i)]
         set yStart [expr $yLeftStart + $k*($yRightStart-$yLeftStart)/$numXEle($i)]
         set xEnd [expr $xLeftEnd + $k*($xRightEnd-$xLeftEnd)/$numXEle($i)]
         set yEnd [expr $yLeftEnd + $k*($yRightEnd-$yLeftEnd)/$numXEle($i)]
         set xIncr [expr ($xEnd - $xStart)/ $numYEle($j)]
         set yIncr [expr ($yEnd - $yStart)/ $numYEle($j)]
         set matNum $diag($i,$j,1)
     
         for {set m 0} {$m < $numYEle($j)} {incr m 1} {
            set xDim [expr $xStart + $m*$xIncr]
            set yDim [expr $yStart + $m*$yIncr]

            if { ($j>0) && ($diag($i,[expr $j-1],0)==1) && ($k==$numXEle($i)) && ($m==0) } {
               # add equalDOF for the block below
               set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
               lappend rnode $nodeNum
               set matNum $diag($i,[expr $j-1],2)
               lappend mat $matNum;

               set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
               equalDOF [expr $nodeNum-1] $nodeNum 1 2
               lappend rnode $nodeNum
               set matNum $diag($i,$j,1)
               lappend mat $matNum;
            } else {
               set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
               lappend rnode $nodeNum
               lappend mat $matNum;
            }
        
            if { $diag($i,$j,0)==1 } {
               set z [expr double($numXEle($i))/double($numYEle($j))]
               set z1 [expr double($k)/double($z)-1]
               if { ($m>$z1) && ($m<[expr $z1+1]) } {
                  set xDim [expr $xDim+(double($k)/double($z)-$m)*$xIncr]
                  set yDim [expr $yDim+(double($k)/double($z)-$m)*$yIncr]
                  set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
                  lappend rnode $nodeNum
                  set matNum $diag($i,$j,2)
                  lappend mat $matNum;
               } elseif {$m==[expr $z1+1]} {
                  set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
                  equalDOF [expr $nodeNum-1] $nodeNum 1 2
                  lappend rnode $nodeNum
                  set matNum $diag($i,$j,2)
                  set mat [lreplace $mat end end $matNum]
                  lappend mat $matNum;
               }
            } 
 
            if { $diag($i,$j,0)==2 } {
               set z [expr double($numXEle($i))/double($numYEle($j))]
               set z1 [expr double($k)/double($z)]
               set z2 [expr $numYEle($j)-$m] 
               if { ($z2>$z1) && ($z2<[expr $z1+1]) } {
                  set xDim [expr $xDim+$xIncr-(double($k)/double($z)-$z2+1)*$xIncr]
                  set yDim [expr $yDim+$yIncr-(double($k)/double($z)-$z2+1)*$yIncr]
                  set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
                  lappend rnode $nodeNum
                  set matNum $diag($i,$j,2)
                  lappend mat $matNum;             
               } elseif { $z2==$z1 } {
                  set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
                  equalDOF [expr $nodeNum-1] $nodeNum 1 2
                  lappend rnode $nodeNum
                  set matNum $diag($i,$j,2)
                  lappend mat $matNum;
               }
            } 
         }
       
         # last block
         if { $j == [expr $numYBlk-1] } {
            set xDim $xEnd
            set yDim $yEnd
            set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
            lappend rnode $nodeNum
            lappend mat $matNum;

            if { ($diag($i,$j,0)==1) && ($k == $numXEle($i)) } {
               set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
               equalDOF [expr $nodeNum-1] $nodeNum 1 2
               lappend rnode $nodeNum
               set matNum $diag($i,$j,2)
               set mat [lreplace $mat end end $matNum]
               lappend mat $matNum;
            }
         }   
      }  
      #puts "$lnode"
      #puts "$rnode"
      #puts "$mat"
      for {set nod 0} {$nod < [expr [llength $lnode]-1]} {incr nod 1} {
         set n1 [lindex $lnode $nod]
         set n2 [lindex $rnode $nod]
         set nodp1 [expr $nod+1]
         set n3 [lindex $rnode $nodp1]
         set n4 [lindex $lnode $nodp1]
         set matOpt [lindex $mat $nod]
         set matDes $eleMat($matOpt)
         set grav $gravY($matOpt)
         incr eleNum 1
                                          #   thick  material   maTag    press   mDensity   gravity 
         element quad $eleNum $n1 $n2 $n3 $n4 $thick1   $matDes   $matOpt  0     $massDen   0 $grav 
         puts $eleOutput "$eleNum $n1 $n2 $n3 $n4 $matOpt"    
      }  
      set lnode $rnode
      set rnode {}
      set mat {};
   }   
}

# abutments generation
# user defined data
set leftParallel 0    ;# leftParallel = 1: more divisions parallel to left slope
                      ;# This is better when the tip is closer to the left 
set numMinorDiv 3     ;# number of minor divisions
set xTip $X(228)    
set yTip [expr $yb+8.2]
set matOpt $matter(6)
set numBaseNode 15
set baseNodeNumList [list 84 96 108 120 132 144 156 168 180 192 204 216 228 242 255]
if {[llength $baseNodeNumList] != $numBaseNode} {
   error "abutment generation: insufficient number of elements in baseNodeNumList"
}
for {set i 0} {$i < $numBaseNode} {incr i 1} {
   set baseNodeNum($i) [lindex $baseNodeNumList $i]
}
unset baseNodeNumList

for {set i 0} {$i < $numBaseNode} {incr i 1} {
   set xBase($i) $X($baseNodeNum($i))
   set yBase($i) $Y($baseNodeNum($i))
}


# find base nodes of minor divisions
for {set i 1} {$i < $numMinorDiv} {incr i 1} {
  set nodeInc [expr $i*$numBaseNode/$numMinorDiv]
  set dualNode($i) $nodeInc
}
if {$leftParallel==1} {
  set dualNode(0) [expr $numBaseNode-1]
} else {
  set dualNode(0) 0
} 
# add dual node at dualNode(0)
set nodeNum [addNode $nodeNum $xBase($dualNode(0)) $yBase($dualNode(0)) $nodeOutput] 
equalDOF $nodeNum $baseNodeNum($dualNode(0)) 1 2


# find left/right slopes 
set kLeft [expr ($yTip - $yBase(0))/($xTip - $xBase(0))]
set kRight [expr ($yTip - $yBase([expr $numBaseNode-1]))/($xTip - $xBase([expr $numBaseNode-1]))]

# generates nodes and elements for the abutments above ground surface
set k 0
set x1 $xBase(0)
for {set i 1} {$i < $numBaseNode} {incr i 1} {
    if {($k < [expr $numMinorDiv-1]) && ($dualNode([expr $k+1]) == $i)} {
      set nodeNum [addNode $nodeNum $xBase($i) $yBase($i) $nodeOutput] 
      equalDOF $nodeNum $baseNodeNum($dualNode([expr $k+1])) 1 2
    }

    for {set j 0} {$j <= $k} {incr j 1} {   
      set x1 $xBase($dualNode([expr $k - $j]))
      set xDim [expr ($kLeft*$x1 - $kRight*$xBase($i))/($kLeft-$kRight)]
      set yDim [expr $kLeft*($xDim-$x1) + $yBase($i)] 
      set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 

      if {$j == 0} {
        set n1 $baseNodeNum([expr $i-1])
        if {$dualNode($k) == [expr $i-1]} {
          set n2 [expr $nodeNum - $k - 1]
          set n3 $baseNodeNum($i)
          set n4 $nodeNum 
        } elseif {($k < [expr $numMinorDiv-1]) && ($dualNode([expr $k+1]) == $i)} {
          set n2 $baseNodeNum($i)
          set n3 $nodeNum 
          set n4 [expr $nodeNum - $k - 2] 
        } else {
          set n2 $baseNodeNum($i)
          set n3 $nodeNum 
          set n4 [expr $nodeNum - $k - 1] 
        }
      } elseif {($k < [expr $numMinorDiv-1]) && ($dualNode([expr $k+1]) == $i)} {
        set n1 [expr $nodeNum - $k - 3] 
        set n2 [expr $nodeNum - 1] 
        set n3 $nodeNum
        set n4 [expr $nodeNum - $k - 2] 
      } else {
        set n1 [expr $nodeNum - $k - 2] 
        set n2 [expr $nodeNum - 1] 
        set n3 $nodeNum
        set n4 [expr $nodeNum - $k - 1] 
      }
  
      set matDes $eleMat($matOpt)
      set grav $gravY($matOpt)
      incr eleNum 1
                                        #   thick  material   maTag    press   mDensity   gravity 
      element quad $eleNum $n1 $n2 $n3 $n4 $thick2   $matDes   $matOpt    0   $massDen   0 $grav 
      puts $eleOutput "$eleNum $n1 $n2 $n3 $n4 $matOpt"    
    }

    if {($k < [expr $numMinorDiv-1]) && ($dualNode([expr $k+1]) == $i)} {
      incr k 1
    }
}

# right abutment
set leftParallel 1    ;# leftParallel = 1: more divisions parallel to left slope
                      ;# This is better when the tip is closer to the left 
set numMinorDiv 3     ;# number of minor divisions
set matOpt $matter(6)
set xTip $X(909)    
set yTip [expr $yb+9.8]

set numBaseNode 15
set baseNodeNumList [list 883 896 909 921 933 945 957 969 981 993 1005 1017 1029 1041 1053]
if {[llength $baseNodeNumList] != $numBaseNode} {
   error "abutment generation: insufficient number of elements in baseNodeNumList"
}
for {set i 0} {$i < $numBaseNode} {incr i 1} {
   set baseNodeNum($i) [lindex $baseNodeNumList $i]
}
unset baseNodeNumList

for {set i 0} {$i < $numBaseNode} {incr i 1} {
   set xBase($i) $X($baseNodeNum($i))
   set yBase($i) $Y($baseNodeNum($i))
}

# find base nodes at minor divisions
for {set i 0} {$i < [expr $numMinorDiv-1]} {incr i 1} {
  set nodeInc [expr ($i+1)*$numBaseNode/$numMinorDiv]
  set dualNode($i) $nodeInc
}
if {$leftParallel==1} {
  set dualNode([expr $numMinorDiv-1]) [expr $numBaseNode-1]
} else {
  set dualNode(0) 0
} 

 
# find left/right slopes 
set kLeft [expr ($yTip - $yBase(0))/($xTip - $xBase(0))]
set kRight [expr ($yTip - $yBase([expr $numBaseNode-1]))/($xTip - $xBase([expr $numBaseNode-1]))]

# generates nodes and elements for the abutments above ground surface
set k [expr $numMinorDiv-1]
for {set i 0} {$i < $numBaseNode} {incr i 1} {
    for {set j 0} {$j <= $k} {incr j 1} {   
      if {($j == 0) && ($dualNode([expr $numMinorDiv-1 - $k]) == $i)} {
        set nodeNum [addNode $nodeNum $xBase($i) $yBase($i) $nodeOutput] 
        equalDOF $nodeNum $baseNodeNum($dualNode([expr $numMinorDiv-1 - $k])) 1 2
      } else {
        set x1 $xBase($dualNode([expr $j + $numMinorDiv-1 - $k]))
        set xDim [expr ($kLeft*$xBase($i) - $kRight*$x1)/($kLeft-$kRight)]
        set yDim [expr $kLeft*($xDim-$xBase($i)) + $yBase($i)] 
        set nodeNum [addNode $nodeNum $xDim $yDim $nodeOutput] 
      }

      if {$i > 0} {
        if {$j == 0} {
          set n1 $baseNodeNum([expr $i-1])
          set n2 $baseNodeNum($i)
          set n3 $nodeNum 
          set n4 [expr $nodeNum - $k - 1] 
        } else {
          set n1 [expr $nodeNum - $k - 2] 
          set n2 [expr $nodeNum - 1] 
          set n3 $nodeNum
          set n4 [expr $nodeNum - $k - 1] 
        }
  
        set matDes $eleMat($matOpt)
        set grav $gravY($matOpt)
        incr eleNum 1
                                          #   thick  material   maTag    press   mDensity   gravity 
        element quad $eleNum $n1 $n2 $n3 $n4 $thick2   $matDes   $matOpt    0     $massDen   0 $grav 
        puts $eleOutput "$eleNum $n1 $n2 $n3 $n4 $matOpt"    
      }
    }
    if {$dualNode([expr $numMinorDiv-1 -$k]) == $i} {
      incr k -1
    }
}

close $nodeOutput 
close $eleOutput 

#boundary conditions
set k 0
for {set i 1} {$i <= $nodeNum} {incr i 1} {   
   if {$X($i)==$xRightMost} { 
     incr k 1
     equalDOF $i $k  1
   }
   if {$Y($i)==0.0} { fix $i 1 1 }
}

unset lnode rnode mat matter diag XList YList numXEle numYEle baseNodeNum

#############################################################
# GRAVITY APPLICATION (static, elastic behavior)

# create the SOE, ConstraintHandler, Integrator, Algorithm and Numberer

test NormUnbalance 1.0e-4 10 0
system ProfileSPD
algorithm Linear 
constraints Penalty 1.e10 1.e10
integrator LoadControl 1 1 1 1
numberer RCM

# create the Analysis
analysis Static 

analyze 2

# switch material stage from elastic (gravity) to plastic
    for {set i 1} {$i <=11 } {incr i 1} {
       updateMaterialStage -material $i -stage 1
    }

###############################################################
# BUILD PILE MODEL

# Define model builder for pile
model basic -ndm 2 -ndf 3

geomTransf Linear 1
set nodeOutput [open pileNode w]  ;# node number and coordinates output
set eleOutput [open pileEle w]    ;# element number and nodes output

# define the nodes
set capElev [expr $yb+1] 
set numPile 8
set capElevaList [list $capElev $capElev $capElev $capElev $capElev $capElev $capElev $capElev]
set pileDepthList [list -16.8 -13.7 -16.8 -12.8 -16.8 -16.8 -12.2 -18.3]  ;# in meters
array set pnod {
  0  308
  1  382
  2  454
  3  526
  4  598
  5  670
  6  750
  7  830
}

array set pileI {
   1   2.7
   2   12.
   3   12.
   4   12.
   5   12.
   6   12.
   7   2.7
   8   2.7
}
set pileE 32000000.
array set pileA {
   1   2.02
   2   7.39
   3   7.39
   4   7.39
   5   7.39
   6   7.39
   7   2.02
   8   2.02
}

for {set i 0} {$i < $numPile} {incr i 1} {
  # find surface node and coordinates
  set capEleva [lindex $capElevaList $i] 
  set x $X($pnod($i))
  set yTip [expr $capEleva + [lindex $pileDepthList $i]]

  incr nodeNum 1
  node $nodeNum $x $capEleva 
  set X($nodeNum) $x ; set Y($nodeNum) $capEleva
  puts $nodeOutput "$nodeNum $x $capEleva"
  set capNode($i) $nodeNum
  set y $capEleva 
  set k 0
  while {$y > $yTip} {
    set soilNode [expr $pnod($i)-$k]
    set y $Y($soilNode)
    incr k 1
    if {$y==$Y([expr $soilNode+1])} {  # skip equalDOF soil node
      set soilNode [expr $pnod($i)-$k]
      set y $Y($soilNode)
      incr k 1
    }
    incr nodeNum 1
    node $nodeNum $x $y
    set X($nodeNum) $x; set Y($nodeNum) $y 
    puts $nodeOutput "$nodeNum $x $y"
    equalDOF  $nodeNum   $soilNode      1  2
    incr eleNum 1
      #                          Tag   Inode   Jnode                A   E     I  tranfTag
    element elasticBeamColumn  $eleNum $nodeNum  [expr $nodeNum-1] $pileA([expr $i+1])  $pileE  $pileI([expr $i+1])  1
    puts $eleOutput "$eleNum $nodeNum  [expr $nodeNum-1] 1"
  } 
}

close $nodeOutput
close $eleOutput

#-------------------------------------UCLA---------------------------------
###############################################################
# Build bridge Model
#------------------------------------------
# Start of Model Generation for the bridge
#------------------------------------------
#loadConst -time 0.0
setTime 0.0

model basic -ndm 2 -ndf 2
set nodeTemp $nodeNumMaxUCSD
incr nodeTemp
set nodeLeft $nodeTemp
incr nodeTemp
set nodeRight $nodeTemp
incr nodeTemp
node $nodeLeft    $X($abutmentLeft)  $Y($abutmentLeft)
set X($nodeLeft) $X($abutmentLeft)
set Y($nodeLeft) $Y($abutmentLeft)
node $nodeRight   $X($abutmentRight) $Y($abutmentRight)
set X($nodeRight) $X($abutmentRight)
set Y($nodeRight) $Y($abutmentRight)


# Create ModelBuilder (with two-dimensions and 3 DOF/node)

model basic -ndm 2 -ndf 3

set beam1LeftNode $nodeTemp
incr nodeTemp
set beam1RightNode $nodeTemp
incr nodeTemp
set col1TopNode $nodeTemp
incr nodeTemp
set beam2RightNode $nodeTemp
incr nodeTemp
set col2TopNode $nodeTemp
incr nodeTemp
set beam3RightNode $nodeTemp
incr nodeTemp
set col3TopNode $nodeTemp
incr nodeTemp
set beam4LeftNode $nodeTemp
incr nodeTemp
set beam4RightNode $nodeTemp
incr nodeTemp
set col4TopNode $nodeTemp
incr nodeTemp
set beam5RightNode $nodeTemp
incr nodeTemp
set col5TopNode $nodeTemp
incr nodeTemp
set beam6RightNode $nodeTemp
incr nodeTemp
set col6TopNode $nodeTemp
incr nodeTemp
set beam7LeftNode $nodeTemp
incr nodeTemp
set beam7RightNode $nodeTemp
incr nodeTemp
set col7TopNode $nodeTemp
incr nodeTemp
set beam8RightNode $nodeTemp
incr nodeTemp
set col8TopNode $nodeTemp
incr nodeTemp
set beam9RightNode $nodeTemp


unset nodeTemp


# Create nodes
#      tag                  x                  y 
node $beam1LeftNode  $X($abutmentLeft)  $Y($abutmentLeft)
set X($beam1LeftNode) $X($abutmentLeft)
set Y($beam1LeftNode) $Y($abutmentLeft)

node $beam1RightNode $X($col1BaseNode)  [expr $Y($col1BaseNode)+$colLength1]
set X($beam1RightNode) $X($col1BaseNode)
set Y($beam1RightNode) [expr $Y($col1BaseNode)+$colLength1]
node $col1TopNode    $X($col1BaseNode)  [expr $Y($col1BaseNode)+$colLength1]
set X($col1TopNode) $X($col1BaseNode)
set Y($col1TopNode) [expr $Y($col1BaseNode)+$colLength1]

node $beam2RightNode $X($col2BaseNode)  [expr $Y($col2BaseNode)+$colLength2]
set X($beam2RightNode) $X($col2BaseNode)
set Y($beam2RightNode) [expr $Y($col2BaseNode)+$colLength2]
node $col2TopNode    $X($col2BaseNode)  [expr $Y($col2BaseNode)+$colLength2]
set X($col2TopNode) $X($col2BaseNode)
set Y($col2TopNode) [expr $Y($col2BaseNode)+$colLength2]

node $beam3RightNode $X($col3BaseNode)  [expr $Y($col3BaseNode)+$colLength3]
set X($beam3RightNode) $X($col3BaseNode)
set Y($beam3RightNode) [expr $Y($col3BaseNode)+$colLength3]
node $col3TopNode    $X($col3BaseNode)  [expr $Y($col3BaseNode)+$colLength3]
set X($col3TopNode) $X($col3BaseNode)
set Y($col3TopNode) [expr $Y($col3BaseNode)+$colLength3]
node $beam4LeftNode  $X($col3BaseNode)  [expr $Y($col3BaseNode)+$colLength3]
set X($beam4LeftNode) $X($col3BaseNode)
set Y($beam4LeftNode) [expr $Y($col3BaseNode)+$colLength3]

node $beam4RightNode $X($col4BaseNode)  [expr $Y($col4BaseNode)+$colLength4]
set X($beam4RightNode) $X($col4BaseNode)
set Y($beam4RightNode) [expr $Y($col4BaseNode)+$colLength4]
node $col4TopNode    $X($col4BaseNode)  [expr $Y($col4BaseNode)+$colLength4]
set X($col4TopNode) $X($col4BaseNode)
set Y($col4TopNode) [expr $Y($col4BaseNode)+$colLength4]

node $beam5RightNode $X($col5BaseNode)  [expr $Y($col5BaseNode)+$colLength5]
set X($beam5RightNode) $X($col5BaseNode)
set Y($beam5RightNode) [expr $Y($col5BaseNode)+$colLength5]
node $col5TopNode    $X($col5BaseNode)  [expr $Y($col5BaseNode)+$colLength5]
set X($col5TopNode) $X($col5BaseNode)
set Y($col5TopNode) [expr $Y($col5BaseNode)+$colLength5]

node $beam6RightNode $X($col6BaseNode)  [expr $Y($col6BaseNode)+$colLength6]
set X($beam6RightNode) $X($col6BaseNode)
set Y($beam6RightNode) [expr $Y($col6BaseNode)+$colLength6]
node $col6TopNode    $X($col6BaseNode)  [expr $Y($col6BaseNode)+$colLength6]
set X($col6TopNode) $X($col6BaseNode)
set Y($col6TopNode) [expr $Y($col6BaseNode)+$colLength6]
node $beam7LeftNode  $X($col6BaseNode)  [expr $Y($col6BaseNode)+$colLength6]
set X($beam7LeftNode) $X($col6BaseNode)
set Y($beam7LeftNode) [expr $Y($col6BaseNode)+$colLength6]

node $beam7RightNode $X($col7BaseNode)  [expr $Y($col7BaseNode)+$colLength7]
set X($beam7RightNode) $X($col7BaseNode)
set Y($beam7RightNode) [expr $Y($col7BaseNode)+$colLength7]
node $col7TopNode    $X($col7BaseNode)  [expr $Y($col7BaseNode)+$colLength7]
set X($col7TopNode) $X($col7BaseNode)
set Y($col7TopNode) [expr $Y($col7BaseNode)+$colLength7]

node $beam8RightNode $X($col8BaseNode)  [expr $Y($col8BaseNode)+$colLength8]
set X($beam8RightNode) $X($col8BaseNode)
set Y($beam8RightNode) [expr $Y($col8BaseNode)+$colLength8]
node $col8TopNode    $X($col8BaseNode)  [expr $Y($col8BaseNode)+$colLength8]
set X($col8TopNode) $X($col8BaseNode)
set Y($col8TopNode) [expr $Y($col8BaseNode)+$colLength8]

node $beam9RightNode $X($abutmentRight) $Y($abutmentRight)
set X($beam9RightNode) $X($abutmentRight)
set Y($beam9RightNode) $Y($abutmentRight)


set nodeOutput [open bridgeNode w]
for {set temp [expr $nodeNumMaxUCSD+1]} {$temp <= $beam9RightNode} {incr temp} {
  puts $nodeOutput "$temp $X($temp) $Y($temp)"
}
close $nodeOutput

set filetemp [open allNodes w]
for { set temp 1 } { $temp <= $beam9RightNode } {incr temp} {
  puts $filetemp "$temp $X($temp) $Y($temp)"
}
close $filetemp 


#              ndR             ndC      dof1? dof2? dof3?
equalDOF  $abutmentLeft   $nodeLeft         2
equalDOF  $nodeLeft       $beam1LeftNode    1    2
equalDOF  $col1TopNode    $beam1RightNode   1    2
equalDOF  $col2TopNode    $beam2RightNode   1    2
equalDOF  $col3TopNode    $beam3RightNode   2
equalDOF  $col3TopNode    $beam4LeftNode    1    2
equalDOF  $col4TopNode    $beam4RightNode   1    2
equalDOF  $col5TopNode    $beam5RightNode   1    2
equalDOF  $col6TopNode    $beam6RightNode   2
equalDOF  $col6TopNode    $beam7LeftNode    1    2
equalDOF  $col7TopNode    $beam7RightNode   1    2
equalDOF  $col8TopNode    $beam8RightNode   1    2
equalDOF  $beam9RightNode $nodeRight        1    2
equalDOF  $nodeRight      $abutmentRight    2


set bridgeMatNum 10
set temp 1
for {set matTemp [expr $matNumMaxUCSD+1]} {$matTemp <= [expr $matNumMaxUCSD+$bridgeMatNum]} {incr matTemp} {
    set bridgeMat$temp $matTemp
    incr temp
}

unset temp matTemp

#Define materials for nonlinear columns
#Core concrete (confined)
#CONCRETE                      tag        f'c     ec0    f'cu     ecu
uniaxialMaterial Concrete01 $bridgeMat1 $fcCore -0.004 $fcuCore -0.014

#Cover concrete (unconfined)
#CONCRETE                      tag         f'c     ec0     f'cu     ecu
uniaxialMaterial Concrete01 $bridgeMat2 $fcCover -0.002 $fcuCover -0.008 

#Reinfocing steel
#STEEL                        tag     fy E0 b
uniaxialMaterial Steel01 $bridgeMat3 $fy $E 0.008 

#Gap                               tag                          E                                     fy                         gap
uniaxialMaterial ElasticPPGap $bridgeMat4   [expr 830530.0*$kips2KN/$in2m]  [expr -386.9*$kips2KN/($in2m*$in2m)]  [expr -0.5*$in2m] ;  # gap
uniaxialMaterial ElasticPPGap $bridgeMat5   [expr 830530.0*$kips2KN/$in2m]  [expr  386.9*$kips2KN/($in2m*$in2m)]  [expr  0.5*$in2m] ;  # hook
uniaxialMaterial Parallel $bridgeMat6 $bridgeMat4 $bridgeMat5

uniaxialMaterial ElasticPPGap $bridgeMat7   [expr 7500.0*$kips2KN/$in2m]  [expr -26760.0*$kips2KN/($in2m*$in2m)] [expr -2.0*$in2m]  ;  # gap

uniaxialMaterial ElasticPPGap $bridgeMat8  [expr 1418700.0*$kips2KN/$in2m] [expr -571.2*$kips2KN/($in2m*$in2m)]  [expr -0.5*$in2m]  ; # gap
uniaxialMaterial ElasticPPGap $bridgeMat9  [expr 1418700.0*$kips2KN/$in2m] [expr  494.4*$kips2KN/($in2m*$in2m)]  [expr  0.5*$in2m]  ; # hook
uniaxialMaterial Parallel $bridgeMat10 $bridgeMat8 $bridgeMat9

#--------------------------------------------------------

#Define cross-section for nonlinear columns
#------------------------------------------------
set columnSectionY1 [expr $pierDepth/2.0-$cover]
set columnSectionZ1 [expr $pierWidth/2.0]

section Fiber 1 {
  
        # Create the concrete core fibers
        patch quad $bridgeMat1 10 1 [expr -$columnSectionY1] [expr -$columnSectionZ1] 0.0 [expr -$columnSectionZ1-$triDepth] 0.0 [expr $columnSectionZ1+$triDepth] [expr -$columnSectionY1] $columnSectionZ1
        patch quad $bridgeMat1 1 10 $columnSectionY1 [expr -$columnSectionZ1] $columnSectionY1 $columnSectionZ1 0.0 [expr $columnSectionZ1+$triDepth] 0.0 [expr -$columnSectionZ1-$triDepth]

        #Create the concrete cover fibers
        patch rect $bridgeMat2 1 1 [expr -$columnSectionY1-$cover] [expr -$columnSectionZ1] [expr -$columnSectionY1] $columnSectionZ1
        patch rect $bridgeMat2 1 1 $columnSectionY1 [expr -$columnSectionZ1] [expr $columnSectionY1+$cover] $columnSectionZ1
        patch quad $bridgeMat2 1 5 0.0 [expr -$columnSectionZ1-$triDepth-$cover] 0.0 [expr -$columnSectionZ1-$triDepth] [expr -$columnSectionY1] [expr -$columnSectionZ1] [expr -$columnSectionY1-$cover] [expr -$columnSectionZ1] 
        patch quad $bridgeMat2 1 5 0.0 [expr -$columnSectionZ1-$triDepth] 0.0 [expr -$columnSectionZ1-$triDepth-$cover] [expr $columnSectionY1+$cover] [expr -$columnSectionZ1] $columnSectionY1 [expr -$columnSectionZ1]
        patch quad $bridgeMat2 1 5 [expr -$columnSectionY1-$cover] $columnSectionZ1 [expr -$columnSectionY1] $columnSectionZ1 0.0 [expr $columnSectionZ1+$triDepth] 0.0 [expr $columnSectionZ1+$triDepth+$cover]
        patch quad $bridgeMat2 1 5 $columnSectionY1 $columnSectionZ1 [expr $columnSectionY1+$cover] $columnSectionZ1 0.0 [expr $columnSectionZ1+$triDepth+$cover] 0.0 [expr $columnSectionZ1+$triDepth]

        #Create the reinforcing fibers
        layer straight $bridgeMat3 1 $As 0.0 [expr -$columnSectionZ1-$triDepth]  0.0 [expr -$columnSectionZ1-$triDepth]
        layer straight $bridgeMat3 2 $As [expr -0.25*$triDepth] [expr -$columnSectionZ1-0.75*$triDepth] [expr 0.25*$triDepth] [expr -$columnSectionZ1-0.75*$triDepth] 
        layer straight $bridgeMat3 2 $As [expr -0.5*$triDepth] [expr -$columnSectionZ1-0.5*$triDepth] [expr 0.5*$triDepth] [expr -$columnSectionZ1-0.5*$triDepth] 
        layer straight $bridgeMat3 2 $As [expr -0.75*$triDepth] [expr -$columnSectionZ1-0.25*$triDepth] [expr 0.75*$triDepth] [expr -$columnSectionZ1-0.25*$triDepth] 
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] [expr -$columnSectionZ1] $columnSectionY1 [expr -$columnSectionZ1]
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] [expr -0.8*$columnSectionZ1] $columnSectionY1 [expr -0.8*$columnSectionZ1]
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] [expr -0.6*$columnSectionZ1] $columnSectionY1 [expr -0.6*$columnSectionZ1]
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] [expr -0.4*$columnSectionZ1] $columnSectionY1 [expr -0.4*$columnSectionZ1]
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] [expr -0.2*$columnSectionZ1] $columnSectionY1 [expr -0.2*$columnSectionZ1]
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] [expr -0.0*$columnSectionZ1] $columnSectionY1 [expr 0.0*$columnSectionZ1]          
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] [expr 0.2*$columnSectionZ1] $columnSectionY1 [expr 0.2*$columnSectionZ1]
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] [expr 0.4*$columnSectionZ1] $columnSectionY1 [expr 0.4*$columnSectionZ1]
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] [expr 0.6*$columnSectionZ1] $columnSectionY1 [expr 0.6*$columnSectionZ1]
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] [expr 0.8*$columnSectionZ1] $columnSectionY1 [expr 0.8*$columnSectionZ1]
        layer straight $bridgeMat3 2 $As [expr -$columnSectionY1] $columnSectionZ1 $columnSectionY1 $columnSectionZ1
        layer straight $bridgeMat3 2 $As [expr -0.75*$triDepth] [expr $columnSectionZ1+0.25*$triDepth] [expr 0.75*$triDepth] [expr $columnSectionZ1+0.25*$triDepth] 
        layer straight $bridgeMat3 2 $As [expr -0.5*$triDepth] [expr $columnSectionZ1+0.5*$triDepth] [expr 0.5*$triDepth] [expr $columnSectionZ1+0.5*$triDepth] 
        layer straight $bridgeMat3 2 $As [expr -0.25*$triDepth] [expr $columnSectionZ1+0.75*$triDepth] [expr 0.25*$triDepth] [expr $columnSectionZ1+0.75*$triDepth] 
        layer straight $bridgeMat3 1 $As 0.0 [expr $columnSectionZ1+$triDepth] 0.0 [expr $columnSectionZ1+$triDepth]
}
#----------------------------------------------------------------------------


#Define column elements
#----------------------------
set eleOutput [open bridgeEle w]
set bridgeColTransf [expr $transfNumMaxUCSD + 1]
set bridgeBeamTransf [expr $transfNumMaxUCSD + 2]
set eleTemp [expr $eleNumMaxUCSD +1]
# Geometry of column elements
#                      tag
geomTransf Linear $bridgeColTransf

# Create the columns using Beam-column elements
#                               tag        ndI          ndJ      nsecs secID transfTag
element nonlinearBeamColumn $eleTemp  $col1BaseNode $col1TopNode  $np     1  $bridgeColTransf
puts $eleOutput "$eleTemp  $col1BaseNode $col1TopNode"
incr eleTemp
element nonlinearBeamColumn $eleTemp  $col2BaseNode $col2TopNode  $np     1  $bridgeColTransf
puts $eleOutput "$eleTemp  $col2BaseNode $col2TopNode"
incr eleTemp
element nonlinearBeamColumn $eleTemp  $col3BaseNode $col3TopNode  $np     1  $bridgeColTransf
puts $eleOutput "$eleTemp  $col3BaseNode $col3TopNode"
incr eleTemp
element nonlinearBeamColumn $eleTemp  $col4BaseNode $col4TopNode  $np     1  $bridgeColTransf
puts $eleOutput "$eleTemp  $col4BaseNode $col4TopNode"
incr eleTemp
element nonlinearBeamColumn $eleTemp  $col5BaseNode $col5TopNode  $np     1  $bridgeColTransf
puts $eleOutput "$eleTemp  $col5BaseNode $col5TopNode"
incr eleTemp
element nonlinearBeamColumn $eleTemp  $col6BaseNode $col6TopNode  $np     1  $bridgeColTransf
puts $eleOutput "$eleTemp  $col6BaseNode $col6TopNode"
incr eleTemp
element nonlinearBeamColumn $eleTemp  $col7BaseNode $col7TopNode  $np     1  $bridgeColTransf
puts $eleOutput "$eleTemp  $col7BaseNode $col7TopNode"
incr eleTemp
element nonlinearBeamColumn $eleTemp  $col8BaseNode $col8TopNode  $np     1  $bridgeColTransf
puts $eleOutput "$eleTemp  $col8BaseNode $col8TopNode"

# Define beam element
#---------------------------------------------
geomTransf Linear $bridgeBeamTransf

incr eleTemp
# Creat the beam element
#                           tag         ndI             ndJ                      A              E         Iz              transfTag
element elasticBeamColumn $eleTemp $beam1LeftNode  $beam1RightNode   [expr 4*$beamSectionArea] $beamE  [expr 4*$beamInertia] $bridgeBeamTransf
puts $eleOutput "$eleTemp $beam1LeftNode  $beam1RightNode"
incr eleTemp
element elasticBeamColumn $eleTemp $beam1RightNode $beam2RightNode   [expr 4*$beamSectionArea] $beamE  [expr 4*$beamInertia] $bridgeBeamTransf
puts $eleOutput "$eleTemp $beam1RightNode $beam2RightNode"
incr eleTemp
element elasticBeamColumn $eleTemp $beam2RightNode $beam3RightNode   [expr 4*$beamSectionArea] $beamE  [expr 4*$beamInertia] $bridgeBeamTransf
puts $eleOutput "$eleTemp $beam2RightNode $beam3RightNode"
incr eleTemp
element elasticBeamColumn $eleTemp $beam4LeftNode  $beam4RightNode   [expr 4*$beamSectionArea] $beamE  [expr 4*$beamInertia] $bridgeBeamTransf
puts $eleOutput "$eleTemp $beam4LeftNode  $beam4RightNode"
incr eleTemp
element elasticBeamColumn $eleTemp $beam4RightNode $beam5RightNode   [expr 4*$beamSectionArea] $beamE  [expr 4*$beamInertia] $bridgeBeamTransf
puts $eleOutput "$eleTemp $beam4RightNode $beam5RightNode"
incr eleTemp
element elasticBeamColumn $eleTemp $beam5RightNode $beam6RightNode   [expr 4*$beamSectionArea] $beamE  [expr 4*$beamInertia] $bridgeBeamTransf
puts $eleOutput "$eleTemp $beam5RightNode $beam6RightNode"
incr eleTemp
element elasticBeamColumn $eleTemp $beam7LeftNode  $beam7RightNode   [expr 4*$beamSectionArea] $beamE  [expr 4*$beamInertia] $bridgeBeamTransf
puts $eleOutput "$eleTemp $beam7LeftNode  $beam7RightNode"
incr eleTemp
element elasticBeamColumn $eleTemp $beam7RightNode $beam8RightNode   [expr 4*$beamSectionArea] $beamE  [expr 4*$beamInertia] $bridgeBeamTransf
puts $eleOutput "$eleTemp $beam7RightNode $beam8RightNode"
incr eleTemp
element elasticBeamColumn $eleTemp $beam8RightNode $beam9RightNode   [expr 4*$beamSectionArea] $beamE  [expr 4*$beamInertia] $bridgeBeamTransf
puts $eleOutput "$eleTemp $beam8RightNode $beam9RightNode"

incr eleTemp
# Define gap element: zeroLength
#       zeroLength    tag       ndI           ndJ         -mat matTag -dir 1
element zeroLength $eleTemp $abutmentLeft   $nodeLeft   -mat $bridgeMat10  -dir 1
puts $eleOutput "$eleTemp $abutmentLeft $nodeLeft"
set abutLGap $eleTemp
incr eleTemp 
element zeroLength $eleTemp $beam3RightNode $col3TopNode    -mat $bridgeMat6   -dir 1
puts $eleOutput "$eleTemp $beam3RightNode $col3TopNode"
set expnLGap $eleTemp
incr eleTemp
element zeroLength $eleTemp $beam6RightNode $col6TopNode    -mat $bridgeMat6   -dir 1
puts $eleOutput "$eleTemp $beam6RightNode $col6TopNode"
set expnRGap $eleTemp
incr eleTemp
element zeroLength $eleTemp $nodeRight      $abutmentRight  -mat $bridgeMat10  -dir 1
puts $eleOutput "$eleTemp $nodeRight $abutmentRight"
set abutRGap $eleTemp
incr eleTemp
element zeroLength $eleTemp $beam3RightNode $beam4LeftNode  -mat $bridgeMat7   -dir 1
puts $eleOutput "$eleTemp $beam3RightNode $beam4LeftNode"
incr eleTemp
element zeroLength $eleTemp $beam6RightNode $beam7LeftNode  -mat $bridgeMat7   -dir 1
puts $eleOutput "$eleTemp $beam6RightNode $beam7LeftNode"

close $eleOutput

set beam1Length [expr $X($col1BaseNode)-$X($abutmentLeft)]
set beam2Length [expr $X($col2BaseNode)-$X($col1BaseNode)]
set beam3Length [expr $X($col3BaseNode)-$X($col2BaseNode)]
set beam4Length [expr $X($col4BaseNode)-$X($col3BaseNode)]
set beam5Length [expr $X($col5BaseNode)-$X($col4BaseNode)]
set beam6Length [expr $X($col6BaseNode)-$X($col5BaseNode)]
set beam7Length [expr $X($col7BaseNode)-$X($col6BaseNode)]
set beam8Length [expr $X($col8BaseNode)-$X($col7BaseNode)]
set beam9Length [expr $X($abutmentRight)-$X($col8BaseNode)]


set gravityForce2 [expr -$g*$wdc*4.0*$beamSectionArea*$beam1Length/2.0]
set gravityForce3 [expr -$g*$wdc*4.0*$beamSectionArea*($beam1Length+$beam2Length)/2.0]
set gravityForce4 [expr -$g*$wdc*$columnSectionArea*0.5*$colLength1]
set gravityForce6 [expr -$g*$wdc*4.0*$beamSectionArea*($beam2Length+$beam3Length)/2.0]
set gravityForce7 [expr -$g*$wdc*$columnSectionArea*0.5*$colLength2]
set gravityForce9 [expr -$g*$wdc*4.0*$beamSectionArea*$beam3Length/2.0]
set gravityForce10 [expr -$g*$wdc*$columnSectionArea*0.5*$colLength3]
set gravityForce12 [expr -$g*$wdc*4.0*$beamSectionArea*$beam4Length/2.0]
set gravityForce13 [expr -$g*$wdc*4.0*$beamSectionArea*($beam4Length+$beam5Length)/2.0]
set gravityForce14 [expr -$g*$wdc*$columnSectionArea*0.5*$colLength4]
set gravityForce16 [expr -$g*$wdc*4.0*$beamSectionArea*($beam5Length+$beam6Length)/2.0]
set gravityForce17 [expr -$g*$wdc*$columnSectionArea*0.5*$colLength5]
set gravityForce19 [expr -$g*$wdc*4.0*$beamSectionArea*$beam6Length/2.0]
set gravityForce20 [expr -$g*$wdc*$columnSectionArea*0.5*$colLength6]
set gravityForce22 [expr -$g*$wdc*4.0*$beamSectionArea*$beam7Length/2.0]
set gravityForce23 [expr -$g*$wdc*4.0*$beamSectionArea*($beam7Length+$beam8Length)/2.0]
set gravityForce24 [expr -$g*$wdc*$columnSectionArea*0.5*$colLength7]
set gravityForce26 [expr -$g*$wdc*4.0*$beamSectionArea*($beam8Length+$beam9Length)/2.0]
set gravityForce27 [expr -$g*$wdc*$columnSectionArea*0.5*$colLength8]
set gravityForce29 [expr -$g*$wdc*4.0*$beamSectionArea*$beam9Length/2.0]

# Define gravity loads
#-----------------------------
pattern Plain 1 "Linear" {

   # Create nodal loads at nodes 
   #       node           FX      FY         MZ
  load $beam1LeftNode   0.0 $gravityForce2  0.0
  load $beam1RightNode  0.0 $gravityForce3  0.0
  load $col1TopNode     0.0 $gravityForce4  0.0
  load $beam2RightNode  0.0 $gravityForce6  0.0
  load $col2TopNode     0.0 $gravityForce7  0.0
  load $beam3RightNode  0.0 $gravityForce9  0.0
  load $col3TopNode     0.0 $gravityForce10 0.0
  load $beam4LeftNode   0.0 $gravityForce12 0.0
  load $beam4RightNode  0.0 $gravityForce13 0.0
  load $col4TopNode     0.0 $gravityForce14 0.0
  load $beam5RightNode  0.0 $gravityForce16 0.0
  load $col5TopNode     0.0 $gravityForce17 0.0
  load $beam6RightNode  0.0 $gravityForce19 0.0
  load $col6TopNode     0.0 $gravityForce20 0.0
  load $beam7LeftNode   0.0 $gravityForce22 0.0
  load $beam7RightNode  0.0 $gravityForce23 0.0
  load $col7TopNode     0.0 $gravityForce24 0.0
  load $beam8RightNode  0.0 $gravityForce26 0.0
  load $col8TopNode     0.0 $gravityForce27 0.0
  load $beam9RightNode  0.0 $gravityForce29 0.0

}

set bridgem2 [expr -$gravityForce2/$g]
set bridgem3 [expr -$gravityForce3/$g]
set bridgem4 [expr -$gravityForce4/$g]
set bridgem6 [expr -$gravityForce6/$g]
set bridgem7 [expr -$gravityForce7/$g]
set bridgem9 [expr -$gravityForce9/$g]
set bridgem10 [expr -$gravityForce10/$g]
set bridgem12 [expr -$gravityForce12/$g]
set bridgem13 [expr -$gravityForce13/$g]
set bridgem14 [expr -$gravityForce14/$g]
set bridgem16 [expr -$gravityForce16/$g]
set bridgem17 [expr -$gravityForce17/$g]
set bridgem19 [expr -$gravityForce19/$g]
set bridgem20 [expr -$gravityForce20/$g]
set bridgem22 [expr -$gravityForce22/$g]
set bridgem23 [expr -$gravityForce23/$g]
set bridgem24 [expr -$gravityForce24/$g]
set bridgem26 [expr -$gravityForce26/$g]
set bridgem27 [expr -$gravityForce27/$g]
set bridgem29 [expr -$gravityForce29/$g]


#          tag            MX        MY      RZ
mass $beam1LeftNode  $bridgem2  $bridgem2  0.0
mass $beam1RightNode $bridgem3  $bridgem3  0.0
mass $col1TopNode    $bridgem4  $bridgem4  0.0
mass $beam2RightNode $bridgem6  $bridgem6  0.0
mass $col2TopNode    $bridgem7  $bridgem7  0.0
mass $beam3RightNode $bridgem9  $bridgem9  0.0
mass $col3TopNode    $bridgem10 $bridgem10 0.0
mass $beam4LeftNode  $bridgem12 $bridgem12 0.0
mass $beam4RightNode $bridgem13 $bridgem13 0.0
mass $col4TopNode    $bridgem14 $bridgem14 0.0
mass $beam5RightNode $bridgem16 $bridgem16 0.0
mass $col5TopNode    $bridgem17 $bridgem17 0.0
mass $beam6RightNode $bridgem19 $bridgem19 0.0
mass $col6TopNode    $bridgem20 $bridgem20 0.0
mass $beam7LeftNode  $bridgem22 $bridgem22 0.0
mass $beam7RightNode $bridgem23 $bridgem23 0.0
mass $col7TopNode    $bridgem24 $bridgem24 0.0
mass $beam8RightNode $bridgem26 $bridgem26 0.0
mass $col8TopNode    $bridgem27 $bridgem27 0.0
mass $beam9RightNode $bridgem29 $bridgem29 0.0

#----------------------------------
# End of the bridge model 
#----------------------------------
#############################################################
# Now apply bridge gravity

wipeAnalysis
constraints Penalty 1.0e10 1.0e10
system ProfileSPD
numberer RCM
test NormDispIncr 5.0e-3 200 1
#test EnergyIncr 1.0e-5 20 1
algorithm ModifiedNewton
integrator LoadControl 0.02 1 0.02 0.02
analysis Static

set startT [clock seconds]
analyze 5
set endT [clock seconds]
puts "Execution time: [expr $endT-$startT] seconds."

puts "bridge gravity analysis completed"

#############################################################
# NOW APPLY LOADING SEQUENCE AND ANALYZE (plastic)

# rezero time
loadConst -time 0.0
#setTime 0.0

wipeAnalysis

#                        patternTag dir   
pattern UniformExcitation    2       1  -accel "Series -factor $accMul -filePath $accNam -dt $accDt"

#recorder for output

recorder Element 482 -time -file stress1.out material 1 stress
recorder Element 482 -time -file strain1.out material 1 strain
recorder Element 484 -time -file stress2.out material 1 stress
recorder Element 484 -time -file strain2.out material 1 strain
recorder Element 486 -time -file stress3.out material 1 stress
recorder Element 486 -time -file strain3.out material 1 strain
recorder Element 488 -time -file stress4.out material 1 stress
recorder Element 488 -time -file strain4.out material 1 strain
recorder Element 490 -time -file stress5.out material 1 stress
recorder Element 490 -time -file strain5.out material 1 strain
recorder Element 492 -time -file stress6.out material 1 stress
recorder Element 492 -time -file strain6.out material 1 strain

recorder Element 1123 -time -file col1.out force
recorder Node col1TopDisp.out disp -time -node $col1TopNode -dof 1 2
recorder Node col1BaseDisp.out disp -time -node $col1BaseNode -dof 1 2
recorder Element 1123 -time -file sn1_strn.out section 1 fiber $columnSectionY1 $columnSectionZ1 strain
recorder Element 1123 -time -file sn1_strs.out section 1 fiber $columnSectionY1 $columnSectionZ1 stress
recorder Element 1123 -time -file cn1_strn.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1123 -time -file cn1_strs.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1123 -time -file cvn1_strn.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 strain
recorder Element 1123 -time -file cvn1_strs.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 stress
recorder Element 1123 -time -file sp1_strn.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 strain
recorder Element 1123 -time -file sp1_strs.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 stress
recorder Element 1123 -time -file cp1_strn.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1123 -time -file cp1_strs.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1123 -time -file cvp1_strn.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 strain
recorder Element 1123 -time -file cvp1_strs.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 stress
recorder Element 1123 -time -file col1SecF.out section 1 force
recorder Element 1123 -time -file col1SecD.out section 1 deformation

recorder Element 1124 -time -file col2.out force
recorder Node col2TopDisp.out disp -time -node $col2TopNode -dof 1 2
recorder Node col2BaseDisp.out disp -time -node $col2BaseNode -dof 1 2
recorder Element 1124 -time -file sn2_strn.out section 1 fiber $columnSectionY1 $columnSectionZ1 strain
recorder Element 1124 -time -file sn2_strs.out section 1 fiber $columnSectionY1 $columnSectionZ1 stress
recorder Element 1124 -time -file cn2_strn.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1124 -time -file cn2_strs.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1124 -time -file cvn2_strn.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 strain
recorder Element 1124 -time -file cvn2_strs.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 stress
recorder Element 1124 -time -file sp2_strn.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 strain
recorder Element 1124 -time -file sp2_strs.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 stress
recorder Element 1124 -time -file cp2_strn.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1124 -time -file cp2_strs.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1124 -time -file cvp2_strn.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 strain
recorder Element 1124 -time -file cvp2_strs.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 stress
recorder Element 1124 -time -file col2SecF.out section 1 force
recorder Element 1124 -time -file col2SecD.out section 1 deformation

recorder Element 1125 -time -file col3.out force
recorder Node col3TopDisp.out disp -time -node $col3TopNode -dof 1 2
recorder Node col3BaseDisp.out disp -time -node $col3BaseNode -dof 1 2
recorder Element 1125 -time -file sn3_strn.out section 1 fiber $columnSectionY1 $columnSectionZ1 strain
recorder Element 1125 -time -file sn3_strs.out section 1 fiber $columnSectionY1 $columnSectionZ1 stress
recorder Element 1125 -time -file cn3_strn.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1125 -time -file cn3_strs.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1125 -time -file cvn3_strn.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 strain
recorder Element 1125 -time -file cvn3_strs.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 stress
recorder Element 1125 -time -file sp3_strn.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 strain
recorder Element 1125 -time -file sp3_strs.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 stress
recorder Element 1125 -time -file cp3_strn.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1125 -time -file cp3_strs.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1125 -time -file cvp3_strn.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 strain
recorder Element 1125 -time -file cvp3_strs.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 stress
recorder Element 1125 -time -file col3SecF.out section 1 force
recorder Element 1125 -time -file col3SecD.out section 1 deformation


recorder Element 1126 -time -file col4.out force
recorder Node col4TopDisp.out disp -time -node $col4TopNode -dof 1 2
recorder Node col4BaseDisp.out disp -time -node $col4BaseNode -dof 1 2
recorder Element 1126 -time -file col4SecD.out section 1 deformation
recorder Element 1126 -time -file col4SecF.out section 1 force 
recorder Element 1126 -time -file sn4_strn.out section 1 fiber $columnSectionY1 $columnSectionZ1 strain
recorder Element 1126 -time -file sn4_strs.out section 1 fiber $columnSectionY1 $columnSectionZ1 stress
recorder Element 1126 -time -file cn4_strn.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1126 -time -file cn4_strs.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1126 -time -file cvn4_strn.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 strain
recorder Element 1126 -time -file cvn4_strs.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 stress
recorder Element 1126 -time -file sp4_strn.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 strain
recorder Element 1126 -time -file sp4_strs.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 stress
recorder Element 1126 -time -file cp4_strn.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1126 -time -file cp4_strs.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1126 -time -file cvp4_strn.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 strain
recorder Element 1126 -time -file cvp4_strs.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 stress

recorder Element 1127 -time -file col5.out force
recorder Node col5TopDisp.out disp -time -node $col5TopNode -dof 1 2
recorder Node col5BaseDisp.out disp -time -node $col5BaseNode -dof 1 2
recorder Element 1127 -time -file sn5_strn.out section 1 fiber $columnSectionY1 $columnSectionZ1 strain
recorder Element 1127 -time -file sn5_strs.out section 1 fiber $columnSectionY1 $columnSectionZ1 stress
recorder Element 1127 -time -file cn5_strn.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1127 -time -file cn5_strs.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1127 -time -file cvn5_strn.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 strain
recorder Element 1127 -time -file cvn5_strs.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 stress
recorder Element 1127 -time -file sp5_strn.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 strain
recorder Element 1127 -time -file sp5_strs.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 stress
recorder Element 1127 -time -file cp5_strn.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1127 -time -file cp5_strs.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1127 -time -file cvp5_strn.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 strain
recorder Element 1127 -time -file cvp5_strs.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 stress
recorder Element 1127 -time -file col5SecF.out section 1 force
recorder Element 1127 -time -file col5SecD.out section 1 deformation

recorder Element 1128 -time -file col6.out force
recorder Node col6TopDisp.out disp -time -node $col6TopNode -dof 1 2
recorder Node col6BaseDisp.out disp -time -node $col6BaseNode -dof 1 2
recorder Element 1128 -time -file sn6_strn.out section 1 fiber $columnSectionY1 $columnSectionZ1 strain
recorder Element 1128 -time -file sn6_strs.out section 1 fiber $columnSectionY1 $columnSectionZ1 stress
recorder Element 1128 -time -file cn6_strn.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1128 -time -file cn6_strs.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1128 -time -file cvn6_strn.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 strain
recorder Element 1128 -time -file cvn6_strs.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 stress
recorder Element 1128 -time -file sp6_strn.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 strain
recorder Element 1128 -time -file sp6_strs.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 stress
recorder Element 1128 -time -file cp6_strn.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1128 -time -file cp6_strs.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1128 -time -file cvp6_strn.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 strain
recorder Element 1128 -time -file cvp6_strs.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 stress
recorder Element 1128 -time -file col6SecF.out section 1 force
recorder Element 1128 -time -file col6SecD.out section 1 deformation

recorder Element 1129 -time -file col7.out force
recorder Node col7TopDisp.out disp -time -node $col7TopNode -dof 1 2
recorder Node col7BaseDisp.out disp -time -node $col7BaseNode -dof 1 2
recorder Element 1129 -time -file col7SecD.out section 1 deformation
recorder Element 1129 -time -file col7SecF.out section 1 force
recorder Element 1129 -time -file sn7_strn.out section 1 fiber $columnSectionY1 $columnSectionZ1 strain
recorder Element 1129 -time -file sn7_strs.out section 1 fiber $columnSectionY1 $columnSectionZ1 stress
recorder Element 1129 -time -file cn7_strn.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1129 -time -file cn7_strs.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1129 -time -file cvn7_strn.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 strain
recorder Element 1129 -time -file cvn7_strs.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 stress
recorder Element 1129 -time -file sp7_strn.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 strain
recorder Element 1129 -time -file sp7_strs.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 stress
recorder Element 1129 -time -file cp7_strn.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1129 -time -file cp7_strs.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1129 -time -file cvp7_strn.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 strain
recorder Element 1129 -time -file cvp7_strs.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 stress

recorder Element 1130 -time -file col8.out force
recorder Node col8TopDisp.out disp -time -node $col8TopNode -dof 1 2
recorder Node col8BaseDisp.out disp -time -node $col8BaseNode -dof 1 2
recorder Element 1130 -time -file sn8_strn.out section 1 fiber $columnSectionY1 $columnSectionZ1 strain
recorder Element 1130 -time -file sn8_strs.out section 1 fiber $columnSectionY1 $columnSectionZ1 stress
recorder Element 1130 -time -file cn8_strn.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1130 -time -file cn8_strs.out section 1 fiber [expr 0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1130 -time -file cvn8_strn.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 strain
recorder Element 1130 -time -file cvn8_strs.out section 1 fiber [expr $columnSectionY1+$cover] $columnSectionZ1 stress
recorder Element 1130 -time -file sp8_strn.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 strain
recorder Element 1130 -time -file sp8_strs.out section 1 fiber [expr -$columnSectionY1] $columnSectionZ1 stress
recorder Element 1130 -time -file cp8_strn.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] strain
recorder Element 1130 -time -file cp8_strs.out section 1 fiber [expr -0.8*$columnSectionY1] [expr 0.1*$columnSectionZ1] stress
recorder Element 1130 -time -file cvp8_strn.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 strain
recorder Element 1130 -time -file cvp8_strs.out section 1 fiber [expr -$columnSectionY1-$cover] $columnSectionZ1 stress
recorder Element 1130 -time -file col8SecF.out section 1 force
recorder Element 1130 -time -file col8SecD.out section 1 deformation


recorder Element $abutLGap -time -file abutLGap.out force
recorder Node abutLeft.out disp -time -node $abutmentLeft -dof 1 2   
recorder Node nodeLeft.out disp -time -node $nodeLeft -dof 1 2

recorder Element $expnLGap -time -file expnLGap.out force
recorder Node beam3RNode.out disp -time -node $beam3RightNode -dof 1 2
recorder Node col3Tnode.out disp -time -node $col3TopNode -dof 1 2 

recorder Element $expnRGap -time -file expnRGap.out force
recorder Node beam6RNode.out disp -time -node $beam6RightNode -dof 1 2
recorder Node col6TNode.out disp -time -node $col6TopNode -dof 1 2

recorder Element $abutRGap -time -file abutRGap.out force
recorder Node nodeRight.out disp -time -node $nodeRight -dof 1 2
recorder Node abutRight.out disp -time -node $abutmentRight -dof 1 2

constraints Penalty 1.0e10 1.0e10
test NormDispIncr 1.e-4 5 1
numberer RCM
#algorithm Newton -initial
algorithm ModifiedNewton -initial
system ProfileSPD
integrator Newmark  $alpha  [expr pow($alpha+0.5, 2)/4]
#analysis Transient
analysis VariableTransient 

initialize

#analyze the structure
set startT [clock seconds]
#analyze $numSteps $deltaT 
analyze $numSteps $deltaT [expr $deltaT/100.] [expr $deltaT/1.] 3
set endT [clock seconds]
puts "Execution time: [expr $endT-$startT] seconds."

wipe  #flush ouput stream
