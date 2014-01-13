# Reliability analysis of steelframe from:
#
if {0} {
@article{Haukaas:2006,
author = {T. Haukaas and M. H. Scott},
title = {Shape Sensitivities in the Reliability Analysis of Nonlinear Frame Structures},
journal = {Computers and Structures},
volume = 84, number = {15--16},
pages = {964--977},
year = 2006
}
}
# OpenSees reliability analysis -- Michael H. Scott (michael.scott@oregonstate.edu)

wipe
wipeReliability

model basic -ndm 2 -ndf 3
reliability

source ~/utility/DiscretizeMember.tcl

set m 1.0
set mm [expr 0.001*$m]
set kN 1.0
set N [expr 0.001*$kN]
set MPa [expr 1000*$kN/($m*$m)]



node 1 [expr  0.0*$m] [expr 0.0*$m]; fix 1 1 1 1
node 2 [expr  5.0*$m] [expr 0.0*$m]; fix 2 1 1 1
node 3 [expr 10.0*$m] [expr 0.0*$m]; fix 3 1 1 1
node 4 [expr 15.0*$m] [expr 0.0*$m]; fix 4 1 1 1

node 5 [expr  0.0*$m] [expr 4.0*$m]
node 6 [expr  5.0*$m] [expr 4.0*$m]
node 7 [expr 10.0*$m] [expr 4.0*$m]
node 8 [expr 15.0*$m] [expr 4.0*$m]

node  9 [expr  0.0*$m] [expr 8.0*$m]
node 10 [expr  5.0*$m] [expr 8.0*$m]
node 11 [expr 10.0*$m] [expr 8.0*$m]
node 12 [expr 15.0*$m] [expr 8.0*$m]

node 13 [expr  0.0*$m] [expr 12.0*$m]
node 14 [expr  5.0*$m] [expr 12.0*$m]
node 15 [expr 10.0*$m] [expr 12.0*$m]
node 16 [expr 15.0*$m] [expr 12.0*$m]

foreach node [getNodeTags] {
   set Y [nodeCoord $node Y]
   set tagY [expr 800+$node]
   randomVariable $tagY normal -mean $Y -stdv [expr 10*$mm]
   parameter $tagY randomVariable $tagY node $node coord 2

   set X [nodeCoord $node X]
   set tagX [expr 900+$node]
   randomVariable $tagX normal -mean $X -stdv [expr 10*$mm+15.0/12*$Y*$mm]
   parameter $tagX randomVariable $tagX node $node coord 1
}


set E [expr 200000.0*$MPa]
set Ecov 0.05

set fy [expr 300.0*$MPa]
set fycov 0.1

set b 0.02
set Hkin [expr $b/(1-$b)*$E]
set Hkincov 0.1

set d [expr 250.0*$mm]
set dcov 0.02

set bf [expr 250.0*$mm]
set bfcov 0.02

set tf [expr 20.0*$mm]
set tfcov 0.02

set tw [expr 20.0*$mm]
set twcov 0.02

uniaxialMaterial Hardening 1 $E $fy 0.0 [expr $b/(1-$b)*$E]
section WFSection2d 1 1 $d $tw $bf $tf 10 2

set nodeTag 20
set eleTag 1

set Nele 4
set Np 4
set integr "Legendre 1 $Np"
geomTransf Linear 1

set members "1 5  5 9  9 13  2 6  6 10  10 14  3 7  7 11  11 15  4 8  8 12  12 16  5 6  9 10  13 14  6 7  10 11  14 15  7 8  11 12  15 16"
set Nmembers 0
foreach {I J} $members {
   incr Nmembers
   set elements($Nmembers) $eleTag
   set eleTag [DiscretizeMember $I $J $Nele dispBeamColumnNL $integr 1 nodeTag $eleTag]
   lappend elements($Nmembers) $eleTag
}


foreach member [array names elements] {

   set e1 [lindex $elements($member) 0]

   set tagE [expr 100+$member]
   randomVariable $tagE lognormal -mean $E -stdv [expr $Ecov*$E]
   parameter $tagE randomVariable $tagE

   set tagfy [expr 200+$member]
   randomVariable $tagfy lognormal -mean $fy -stdv [expr $fycov*$fy]
   parameter $tagfy randomVariable $tagfy

   set tagHkin [expr 300+$member]
   randomVariable $tagHkin lognormal -mean $Hkin -stdv [expr $Hkincov*$Hkin]
   parameter $tagHkin randomVariable $tagHkin

   set tagd [expr 400+$member]
   randomVariable $tagd normal -mean $d -stdv [expr $dcov*$d]
   parameter $tagd randomVariable $tagd

   set tagbf [expr 500+$member]
   randomVariable $tagbf normal -mean $bf -stdv [expr $bfcov*$bf]
   parameter $tagbf randomVariable $tagbf

   set tagtf [expr 600+$member]
   randomVariable $tagtf normal -mean $tf -stdv [expr $tfcov*$tf]
   parameter $tagtf randomVariable $tagtf

   set tagtw [expr 700+$member]
   randomVariable $tagtw normal -mean $tw -stdv [expr $twcov*$tw]
   parameter $tagtw randomVariable $tagtw

   set e2 [lindex $elements($member) 1]

   for {set e $e1} {$e < $e2} {incr e} {
      addToParameter $tagE element $e E
      addToParameter $tagfy element $e fy
      addToParameter $tagHkin element $e Hkin
      addToParameter $tagd element $e d
      addToParameter $tagbf element $e bf
      addToParameter $tagtf element $e tf
      addToParameter $tagtw element $e tw
   }
}

correlateGroup 101 [expr 100+$Nmembers] 0.6
correlateGroup 201 [expr 200+$Nmembers] 0.6
correlateGroup 301 [expr 300+$Nmembers] 0.6

puts "Num parameters: [llength [getParamTags]]"


pattern Plain 1 Linear {
   load 13 [expr 400.0*$kN] 0.0 0.0
   load  9 [expr 267.0*$kN] 0.0 0.0
   load  5 [expr 133.0*$kN] 0.0 0.0
}

pattern Plain 2 Constant {
   foreach node "5 8 9 12 13 16" {
      load $node 0.0 [expr -50.0*$kN] 0.0
   }
   foreach node "6 7 10 11 14 15" {
      load $node 0.0 [expr -100.0*$kN] 0.0
   }
}

set Nsteps 100
integrator LoadControl [expr 1.0/$Nsteps]

analysis Static

sensitivityIntegrator -static
sensitivityAlgorithm -computeAtEachStep






parameter 5000 node 13 disp 1

set Umax [expr 0.03*[nodeCoord 13 Y]]

performanceFunction 1 "$Umax-\$par(5000)"

randomNumberGenerator        CStdLib
probabilityTransformation    Nataf            -print 3
reliabilityConvergenceCheck  Standard         -e1 1.0e-4    -e2 1.0e-4  -print 1
functionEvaluator                Tcl -file "analyze $Nsteps"
gradientEvaluator Implicit
searchDirection              iHLRF
meritFunctionCheck           AdkZhang         -multi 2.0    -add 50    -factor 0.5  
stepSizeRule                 Armijo           -maxNum 10   -base 0.5   -initial 0.3 5
stepSizeRule Fixed -stepSize 1.0
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 20;# -printDesignPointX designPointX.out


runFORMAnalysis steelFrame.out

