# Reliability analysis of industrial frame building from:
#
#@article{Buonopane:2008,
#author = {Stephen G. Buonopane},
#title = {Strength and Reliability of Steel Frames with Random Properties},
#publisher = {ASCE},
#year = {2008},
#journal = {Journal of Structural Engineering},
#volume = {134},
#number = {2},
#pages = {337-344},
#keywords = {Structural reliability; Monte Carlo method; Steel frames; Structural stability; Residual stress}
#}
#
# OpenSees reliability analysis -- Michael H. Scott (michael.scott@oregonstate.edu)

wipe

# Units: kN, m
set in 0.0254
set mm 0.001

model basic -ndm 2 -ndf 3

reliability

# Bay width
set H [expr 10+2.0/3]

# Column height
set L 5.49

# Initial out-of-plumbness
set Delta [expr 1*$L/500]

# Material properties
set E [expr 2.0e8]
set fy [expr 3.45e5]
set alpha 0.00

uniaxialMaterial Elastic 1 $E
uniaxialMaterial Hardening 2 $E $fy 0.0 [expr $alpha*$E/(1-$alpha)]

set matTag 2

# W10x49 Interior columns
section WFSection2d 1 $matTag [expr 9.98*$in] [expr 0.34*$in] [expr 10.0*$in] [expr 0.56*$in] 20 2

# W27x84 Girder
section WFSection2d 2 $matTag [expr 26.71*$in] [expr 0.46*$in] [expr 9.96*$in] [expr 0.64*$in] 20 2

# W12x72 Exterior columns
section WFSection2d 3 $matTag [expr 12.25*$in] [expr 0.43*$in] [expr 12.04*$in] [expr 0.67*$in] 20 2

geomTransf Corotational 1
geomTransf Corotational 2

# Plastic hinge length assumed for interior columns
set lp [expr 1.5*9.98*$in]

set integr "HingeRadauTwo 1 $lp 1 $lp 1"
set Np 6

# PDF assigned to X-coordinates
set Xdist normal; set sigmaX [expr $Delta]; set Xpdf N
#set Xdist uniform; set sigmaX [expr $Delta/sqrt(3.0)]; set Xpdf Uni


node 101 [expr -0.5*$H] 0.0; fix 101 1 1 1
node 102 [expr  0.5*$H] 0.0; fix 102 1 1 1


node 201 [expr -0.5*$H+$Delta] $L

# Make coordinates of node 201 random
randomVariable 201 $Xdist -mean [expr -0.5*$H] -stdv $sigmaX
parameter 201 randomVariable 201 node 201 coord 1
randomVariable 1201 $Xdist -mean $L -stdv $sigmaX
parameter 1201 randomVariable 1201 node 201 coord 2

node 202 [expr  0.5*$H+$Delta] $L

# Make coordinates of node 202 random
randomVariable 202 $Xdist -mean [expr 0.5*$H] -stdv $sigmaX
parameter 202 randomVariable 202 node 202 coord 1
randomVariable 1202 $Xdist -mean $L -stdv $sigmaX
parameter 1202 randomVariable 1202 node 202 coord 2

# Interior girder and column elements
element forceBeamColumn 1 101 201 1 $integr
element forceBeamColumn 2 102 202 1 $integr
element dispBeamColumn     30 201 202 2 Lobatto 2 5



# Exterior column nodes
node 103 [expr -1.5*$H] 0.0; fix 103 1 1 0
node 203 [expr -1.5*$H+$Delta]  $L

node 303 [expr -1.5*$H+$Delta]  $L; equalDOF 203 303 1 2

randomVariable 303 $Xdist -mean [expr -1.5*$H] -stdv $sigmaX
parameter 303 randomVariable 303 node 303 coord 1

randomVariable 1303 $Xdist -mean $L -stdv $sigmaX
parameter 1303 randomVariable 1303 node 303 coord 2

element forceBeamColumn 11 103 303 1 Legendre 3 3
element forceBeamColumn 31 203 201 2 Lobatto 2 5



# Exterior column nodes
node 104 [expr 1.5*$H] 0.0; fix 104 1 1 0
node 204 [expr 1.5*$H+$Delta]  $L

node 304 [expr 1.5*$H+$Delta]  $L; equalDOF 204 304 1 2

randomVariable 304 $Xdist -mean [expr 1.5*$H] -stdv $sigmaX
parameter 304 randomVariable 304 node 304 coord 1

randomVariable 1304 $Xdist -mean $L -stdv $sigmaX
parameter 1304 randomVariable 1304 node 304 coord 2

element forceBeamColumn 21 104 304 1 Legendre 3 3
element forceBeamColumn 41 202 204 2 Lobatto 2 5




# Create group of exterior columns
for {set i 2} {$i <= 5} {incr i} {
    element forceBeamColumn [expr 10+$i] 103 303 1 Legendre 3 3
    element forceBeamColumn [expr 20+$i] 104 304 1 Legendre 3 3
}


# Nominal dead and live distributed loads
set DL 40.9
set LL 20.4

set wn [expr $DL + $LL]
set w [expr 1.2*$DL + 1.6*$LL]

# Nominal dead and live exterior column loads
set DL 1744.0
set LL 872.0

set Pn [expr $DL + $LL]
set P [expr 1.2*$DL + 1.6*$LL]

pattern Plain 1 Linear {

    # Distributed load on girders
    foreach ele [getEleTags] {
	if {$ele >= 30} {
	    eleLoad -ele $ele -type beamUniform -$w
	}
    }

    # Point load on exterior columns
    load 203 0.0 -$P 0.0
    load 204 0.0 -$P 0.0

}

test NormDispIncr 1.0e-8 10 0
algorithm Newton

set Ufinal [expr 80.0*$mm]
set Nsteps 100
set dU [expr $Ufinal/$Nsteps]

integrator DisplacementControl 203 1 $dU

system UmfPack
analysis Static

recorder display "Industrial Frame" 10 10 600 600 -wipe
vup 0 1 0
prp 0 0 200
display 1 0 50






# Assign random variables to interior columns
foreach ele "1 2" {

    # Yield strength
    randomVariable [expr 20+$ele] lognormal -mean [expr 1.*$fy] -stdv [expr 0.06*$fy]
    parameter [expr 20+$ele] randomVariable [expr 20+$ele] element $ele fy
    
    # Elastic modulus
    randomVariable [expr 30+$ele] lognormal -mean [expr 1.*$E] -stdv [expr 0.034*$E]
    parameter [expr 30+$ele] randomVariable [expr 30+$ele] element $ele E

    # Plastic hinge lengths
    randomVariable [expr 500+$ele] normal -mean [expr 1.*$lp] -stdv [expr 0.1*$lp]
    parameter [expr 500+$ele] randomVariable [expr 500+$ele] element $ele lpI

    randomVariable [expr 550+$ele] normal -mean [expr 1.*$lp] -stdv [expr 0.1*$lp]
    parameter [expr 550+$ele] randomVariable [expr 550+$ele] element $ele lpJ

    # Cross-section depth
    randomVariable [expr 600+$ele] normal -mean [expr 9.98*$in] -stdv [expr 0.02*9.98*$in]
    parameter [expr 600+$ele] randomVariable [expr 600+$ele] element $ele d
}

sensitivityIntegrator -static
sensitivityAlgorithm -computeByCommand



set output [open iframe.out w]

puts ""
puts "Frame analysis using mean values of random variables"
puts ""
for {set i 1} {$i <= $Nsteps} {incr i} {
    set ok [analyze 1]
    if {$ok < 0} {
	break
    }

    puts $output "[expr [nodeDisp 203 1]/$mm] [getLoadFactor 1]"
}

close $output

puts "Basic forces of element 1: [basicForce 1]"
puts "Basic forces of element 2: [basicForce 2]"



set Umax [expr 50*$mm]

parameter 1000 node 201 disp 1

performanceFunction 1 "$Umax-\$par(1000)"

reset

set dlam 0.01
set Nsteps [expr int(1.0/$dlam)]

integrator LoadControl $dlam


randomNumberGenerator        CStdLib
probabilityTransformation    Nataf            -print 3
reliabilityConvergenceCheck  Standard         -e1 1.0e-2    -e2 1.0e-2  -print 1
functionEvaluator                Tcl -file "analyze $Nsteps"
gradientEvaluator Implicit
searchDirection              iHLRF
meritFunctionCheck           AdkZhang         -multi 2.0    -add 50    -factor 0.5  
stepSizeRule                 Armijo           -maxNum 10   -base 0.5   -initial 0.3 5
stepSizeRule Fixed -stepSize 1.0
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 20;# -printDesignPointX designPointX.out


runFORMAnalysis industrialFrame.out

puts ""
puts "First order reliability index, beta = [format %.3f $betaFORM(1)]"
puts ""
puts "Importance ranking of random variables at design point"

set remainingRVs [lsort -real [getRVTags]]
set Nrvs [llength $remainingRVs]
for {set i 0} {$i < $Nrvs} {incr i} {
    set largestGamma -1
    foreach rv $remainingRVs {
	set gamma [expr abs($gammaFORM(1,$rv))]
	if {$gamma > $largestGamma} {
	    set largestGamma $gamma
	    set whichRV $rv
	    set whichIndex [lsearch -exact $remainingRVs $rv]
	}
    }

    set xstar $par($whichRV)
    set gamma $gammaFORM(1,$whichRV)
    set mean [getMean $whichRV]
    set stdv [getStdv $whichRV]
    set cov [expr abs($stdv/$mean)]
    puts "RV $whichRV: mean=[format %7.3e $mean], stdev=[format %7.3e $stdv], x*=[format %7.3e $xstar], gamma=[format %7.3e $gamma]"
    set remainingRVs [lreplace $remainingRVs $whichIndex $whichIndex]
}





reset

set Ufinal [expr 80.0*$mm]
set Nsteps 100
set dU [expr $Ufinal/$Nsteps]

integrator DisplacementControl 203 1 $dU

set output [open iframeStar.out w]

puts ""
puts "Frame analysis using design point values of random variables"
puts ""
for {set i 1} {$i <= $Nsteps} {incr i} {
    set ok [analyze 1]
    if {$ok < 0} {
	break
    }
    puts $output "[expr [nodeDisp 203 1]/$mm] [getLoadFactor 1]"
}

close $output

puts "Basic forces of element 1: [basicForce 1]"
puts "Basic forces of element 2: [basicForce 2]"


