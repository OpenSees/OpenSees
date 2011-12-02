# RCFrame3.tcl: R/C two story, two bay frame
# Units: kip, in
# MHS, Sept 1999
#   email: mhscott@ce.berkeley.edu
#
# Elastic analysis for comparison with pushover curve
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
# NOTE: to RUN this example, run the g3 interpreter and 
#       type the command: source RCFrame3.tcl
#
# $Revision: 1.2 $
# $Date: 2002-12-17 02:03:54 $
# $Source: /usr/local/cvs/OpenSees/EXAMPLES/ExampleScripts/RCFrame3.tcl,v $


model BasicBuilder -ndm 2 -ndf 3

#    tag  X   Y
node  1   0   0
node  2   0 180
node  3   0 324
node  4 288   0
node  5 288 180
node  6 288 324
node  7 576   0
node  8 576 180
node  9 576 324

#   node DX DY RZ
fix   1   1  1  1
fix   4   1  1  1
fix   7   1  1  1

# "Cracked" second moments of area
set IcolInt [expr 0.7*1/12*24*pow(27,3)]
set IcolExt [expr 0.7*1/12*24*pow(24,3)]
set Igir    [expr 0.5*1/12*18*pow(24,3)]

# Cross-sectional area of members
set AcolInt [expr 24*27]
set AcolExt [expr 24*24]
set Agir    [expr 18*24]

# Concrete elastic modulus
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

# Gravity load
set P -192
# Lateral load
set H  240

pattern Plain 1 Linear {
   #   node  FX          FY  MZ
   load  2   $H          $P 0.0
   load  5  0.0 [expr $P*2] 0.0
   load  8  0.0          $P 0.0
   load  3   $H [expr $P/2] 0.0
   load  6  0.0          $P 0.0
   load  9  0.0 [expr $P/2] 0.0
}

algorithm Linear
numberer RCM
constraints Plain
integrator LoadControl 1.0
system BandSPD

analysis Static

# Perform the linear analysis
analyze 1

# Get the roof displacement
print node 3


