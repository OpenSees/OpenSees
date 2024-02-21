# RCFrame4.tcl: R/C two story, two bay frame subject to sinusoidal pulse GMotion
# Units: kip, in
# MHS, Sept 1999
#   email: mhscott@ce.berkeley.edu
#
# Linear and Non-linear EQ analysis
#
#    _________________________   _
#   |            |            |  
#   |            |            | 12'
#   |            |            |
#   |____________|____________|  _
#   |            |            |
#   |            |            |
#   |            |            | 15'
#   |            |            |
#   |            |            |  _
#  ===          ===          ===
#   |     24'    |     24'    |
#
#
# NOTE: to RUN this example, run the g3 interpreter and 
#       type the command: source RCFrame5.tcl
#
# $Revision: 1.3 $
# $Date: 2002-12-17 02:03:54 $
# $Source: /usr/local/cvs/OpenSees/EXAMPLES/ExampleScripts/RCFrame5.tcl,v $

model BasicBuilder -ndm 2 -ndf 3

set analysis "LINEAR"
#set analysis "NONLINEAR"

set patternType "UniformExcitation"
#set patternType "MultipleSupport"

#set displayMode "displayON"
set displayMode "displayOFF"

set g 386.4
set P 192
set m [expr $P/$g]
set mR 0
set PI [expr 2.0 * asin(1.0)]
set period 2.5;
set mag 5.0;
set w [expr 2.0 * $PI / $period]

#initial velocity condition for uniform excitation = -udotg
set vel0 [expr -1.0 * $mag * $w]

# ground motion 
set dispSeries "Sine 0 20.0 $period -factor $mag"
set accelSeries "Sine 0 20.0 $period -factor [expr -1.0 * $w * $w * $mag]"

#    tag  X   Y              massX massY rotZ
node  1   0   0
node  2   0 180 -mass           $m          $m      $mR
node  3   0 324 -mass  [expr $m/2] [expr $m/2]      $mR
node  4 288   0
node  5 288 180 -mass  [expr $m/2] [expr $m/2]      $mR
node  6 288 324 -mass           $m          $m      $mR
node  7 576   0
node  8 576 180 -mass           $m          $m      $mR
node  9 576 324 -mass  [expr $m/2] [expr $m/2]      $mR

#   node DX DY RZ
fix   1   1  1  1
fix   4   1  1  1
fix   7   1  1  1

if {$analysis == "NONLINEAR"} {

   #                  tag -f'c  -epsco  -f'cu -epscu
   uniaxialMaterial Concrete01 1 -4.00  -0.002    0.0 -0.006
   uniaxialMaterial Concrete01 2 -5.20  -0.005  -4.70  -0.02

	# Steel model
	#                        tag fy   E     b
	uniaxialMaterial Steel01  3  60 30000 0.02

    # Interior column section
    section fiberSec 1 {
       #           mat nfIJ nfJK   yI  zI    yJ  zJ    yK  zK    yL  zL
       patch quadr  2    1   12 -11.5  10 -11.5 -10  11.5 -10  11.5  10
       patch quadr  1    1   14 -13.5 -10 -13.5 -12  13.5 -12  13.5 -10
       patch quadr  1    1   14 -13.5  12 -13.5  10  13.5  10  13.5  12
       patch quadr  1    1    2 -13.5  10 -13.5 -10 -11.5 -10 -11.5  10
       patch quadr  1    1    2  11.5  10  11.5 -10  13.5 -10  13.5  10
       
       #              mat nBars area    yI zI    yF zF
       layer straight  3    6   1.56 -10.5  9 -10.5 -9
       layer straight  3    6   1.56  10.5  9  10.5 -9
   }

   # Exterior column section
   section fiberSec 2 {
       patch quadr 2 1 10 -10  10 -10 -10  10 -10  10  10
       patch quadr 1 1 12 -12 -10 -12 -12  12 -12  12 -10
       patch quadr 1 1 12 -12  12 -12  10  12  10  12  12
       patch quadr 1 1  2 -12  10 -12 -10 -10 -10 -10  10
       patch quadr 1 1  2  10  10  10 -10  12 -10  12  10
       layer straight 3 6 0.79 -9 9 -9 -9
       layer straight 3 6 0.79  9 9  9 -9
   }

   # Girder section
   section fiberSec 3 {
       patch quadr 1 1 12 -12 9 -12 -9 12 -9 12 9
       layer straight 3 4 1.00 -9 9 -9 -9
       layer straight 3 4 1.00  9 9  9 -9
   }

   set np 4

   # Geometric transformation
   geomTransf Linear 1

   # Define elements
   #                           tag ndI ndJ  nPts  secID transf
   element nonlinearBeamColumn  1   1   2    $np    2      1
   element nonlinearBeamColumn  2   2   3    $np    2      1
   element nonlinearBeamColumn  3   4   5    $np    1      1
   element nonlinearBeamColumn  4   5   6    $np    1      1
   element nonlinearBeamColumn  5   7   8    $np    2      1
   element nonlinearBeamColumn  6   8   9    $np    2      1
   element nonlinearBeamColumn  7   2   5    $np    3      1
   element nonlinearBeamColumn  8   5   8    $np    3      1
   element nonlinearBeamColumn  9   3   6    $np    3      1
   element nonlinearBeamColumn 10   6   9    $np    3      1

   test NormUnbalance 1.0e-06 10 0
   algorithm Newton
#   system SparseGeneral -piv
system Mumps
}

if {$analysis == "LINEAR"} {
   
   # "Cracked" second moments of area
   set IcolInt [expr 0.7*1/12*24*pow(27,3)]
   set IcolExt [expr 0.7*1/12*24*pow(24,3)]
   set Igir    [expr 0.5*1/12*18*pow(24,3)]

   set AcolInt [expr 24*27]
   set AcolExt [expr 24*24]
   set Agir    [expr 18*24]

   set E 4000
   
   # Geometric transformation
   geomTransf Linear 1

   #                         tag ndI ndJ    A     E     I    transf
   element elasticBeamColumn  1   1   2 $AcolExt $E $IcolExt    1
   element elasticBeamColumn  2   2   3 $AcolExt $E $IcolExt    1
   element elasticBeamColumn  3   4   5 $AcolInt $E $IcolInt    1
   element elasticBeamColumn  4   5   6 $AcolInt $E $IcolInt    1
   element elasticBeamColumn  5   7   8 $AcolExt $E $IcolExt    1
   element elasticBeamColumn  6   8   9 $AcolExt $E $IcolExt    1
   element elasticBeamColumn  7   2   5    $Agir $E    $Igir    1
   element elasticBeamColumn  8   5   8    $Agir $E    $Igir    1
   element elasticBeamColumn  9   3   6    $Agir $E    $Igir    1
   element elasticBeamColumn 10   6   9    $Agir $E    $Igir    1

   algorithm Linear
#   system BandSPD
system Mumps
}

if {$patternType == "MultipleSupport"} {
 pattern MultipleSupport  1   {
   groundMotion 1  Series -disp $dispSeries 

   imposedSupportMotion 1 1 1
   imposedSupportMotion 4 1 1
   imposedSupportMotion 7 1 1
  }

  # Create a recorder which writes to Node.out and prints
  # the current time and all dof displacement at node 3 and 2
  recorder Node -file MultiNodeDisp1.out -time -node 3 2 1 -dof 1 2 3 disp
}


if {$patternType == "UniformExcitation"} {
  pattern UniformExcitation  1 1 -accel $accelSeries -vel0 $vel0

  recorder Node -file UniformNodeDisp1.out -time -node 3 2 1 -dof 1 2 3 disp
}

# Source in some g3 commands to display the model
if {$displayMode == "displayON"} {
    # a window to show the displayed shape
    source RCFrameDisplay.tcl 
}

#                   gamma beta
integrator Newmark  0.5  0.25

numberer RCM
constraints Plain

# spit out the 1'st 4 eigen values
eigen 4

# create the analysis object
analysis Transient

#type "Starting Transient Analysis .. hang on"


#        numSteps   dt
analyze    1000    0.01

# spit out the 1'st 4 eigen values
eigen 4




