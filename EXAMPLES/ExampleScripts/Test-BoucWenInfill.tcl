# BoucWenInfill uniaxialMaterial - Example  

# Reference to the experimental campaign:
# ##############################################################################################################
# ##### Cavaleri, L., and F. Di Trapani. 2014. “Cyclic response of masonry infilled    #########################
# ##### RC frames: Experimental results and simplified modeling.” Soil Dyn. Earthquake Eng. 65 (Oct): 224–242. #
# ##### https://doi.org/10.1016/j.soildyn.2014.06.016   ########################################################

# Reference to BoucWenInfill material:
# ##################################################################################################################################
# ##### Sirotti, S., Pelliciari, M., Di Trapani, F., Briseghella, B., Carlo Marano, G., Nuti, C., & Tarantino, A. M. (2021). #######
# ##### Development and validation of new Bouc–Wen data-driven hysteresis model for masonry infilled RC frames. ####################
# ##### Journal of Engineering Mechanics, 147(11), 04021092.  ######################################################################
# ##################################################################################################################################

# Written by: S. Sirotti, Ph.D. Student @ University of Modena and Reggio Emilia, Italy      
# stefano.sirotti@unimore.it       
# Date: February 2022

#Units Measures  kN, mm, sec , kg
wipe

model basic -ndm 2 -ndf 3

#parameters for overall model geometry
set width     4570.0
set height    3125.0

#Create Nodes
#NODE       tag      X      Y      
node         1      0.0    0.0     
node         2    $width   0.0     
node         3      0.0   $height  
node         4    $width  $height  
node         5      0.0   $height

#Constraints at the base of the columns
#fix    tag     Dx      Dy      Rz
fix      1       1       1       0
fix      2       1       1       0
fix      3       0       0       0 
fix      4       0       0       0
fix      5       1       1       1

# Define materials for nonlinear columns
# ------------------------------------------
# CONCRETE                    
# Core concrete (confined) 	   Tag    fpc      epsc0    fpcu     epscu    lambda   ft     Et
uniaxialMaterial Concrete02     1   -0.035    -0.002   -0.028    -0.01     0.10   0.00   0.0
# Cover concrete (unconfined)
uniaxialMaterial Concrete02     2   -0.032    -0.002   -0.025    -0.008    0.10   0.00   0.0

# STEEL
# Reinforcing longitudinal bars 
set fy 0.450;      # Yield stress
set E 210.0;    # Young's modulus
#                        tag   fy  E0    b       R0   cR1   CR1
uniaxialMaterial Steel02  3   $fy  $E  0.1     15  0.925  0.15;

# INFILL                            tag   mass   alpha  beta0   eta0   n    k     xy     deltak  deltaf   psi     Zs    As   epsp   tolerance   MaxNumIter
uniaxialMaterial BoucWenInfill       4    7416   0.003  0.295   0.2   1.5  153  1.8627  0.023    0.0011  0.001   0.17  5.8   0.21   1.0e-05     1000000 

# #########################################################################################################
# set some parameters for COLUMNS sections
set colWidth 350.0
set colDepth 350.0

set cover  30.0
set As    380.0 ;  # area of 1 bar

set y1 [expr $colDepth/2.0]
set z1 [expr $colWidth/2.0]

section Fiber 1 {

    # Create the concrete core fibers
    patch rect 1 30 30 [expr $cover-$y1] [expr $cover-$z1] [expr $y1-$cover] [expr $z1-$cover]

    # Create the concrete cover fibers (top, bottom, left, right)
    patch rect 2 20 20  [expr -$y1] [expr $z1-$cover] $y1 $z1
    patch rect 2 20 20  [expr -$y1] [expr -$z1] $y1 [expr $cover-$z1]
    patch rect 2 20 20  [expr -$y1] [expr $cover-$z1] [expr $cover-$y1] [expr $z1-$cover]
    patch rect 2 20 20  [expr $y1-$cover] [expr $cover-$z1] $y1 [expr $z1-$cover]

    # Create the reinforcing fibers (left, middle, right)
    layer straight 3 3 $As [expr -$y1+$cover] [expr $z1-$cover]  [expr $y1-$cover] [expr $z1-$cover]
	layer straight 3 2 $As [expr -$y1+$cover] 0.0 [expr $y1-$cover] 0.0
    layer straight 3 3 $As [expr -$y1+$cover] [expr -$z1+$cover] [expr $y1-$cover] [expr -$z1+$cover] 
}    

# #################
set Gc 25000000
set C250 10
#______________________
# column torsional stiffness
# Linear elastic torsion for the column

set GJcol [expr $Gc*$C250*$colDepth*pow($colWidth,3)]
set GAcol [expr $Gc*$colWidth*$colDepth*5/6]

uniaxialMaterial Elastic 50 $GJcol
uniaxialMaterial Elastic 51 $GAcol

# Attach torsion to the RC beam section
#                 tag uniTag uniCode secTag
#section Aggregator $secTag $matTag1 $string1 $matTag2 $string2 ....... <-section $sectionTag>
section Aggregator     10    51 Vy      51 Vz    50 T      -section 1 
#______________________
# ##############################################
# set some parameters for BEAM sections
set beaWidth 350.0
set beaDepth 350.0

set cover  30.0
set Ast    153.0;     # area of 1 bar

# some variables derived from the parameters
set yb1 [expr $beaWidth/2.0]
set zb1 [expr $beaDepth/2.0]

section Fiber 4 {

    # Create the concrete core fibers
    patch rect 1 40 40 [expr $cover-$yb1] [expr $cover-$zb1] [expr $yb1-$cover] [expr $zb1-$cover]

    # Create the concrete cover fibers (top, bottom, left, right)
    patch rect 2 20 10  [expr -$yb1] [expr $zb1-$cover] $yb1 $zb1
    patch rect 2 20 10  [expr -$yb1] [expr -$zb1] $yb1 [expr $cover-$zb1]
    patch rect 2 20 10  [expr -$yb1] [expr $cover-$zb1] [expr $cover-$yb1] [expr $zb1-$cover]
    patch rect 2 20 10  [expr $yb1-$cover] [expr $cover-$zb1] $yb1 [expr $zb1-$cover]

    # Create the reinforcing fibers (left, middle, right)
    layer straight 3 4 $Ast [expr -$yb1+$cover] [expr $zb1-$cover] [expr -$yb1+$cover] [expr -$zb1+$cover]
    layer straight 3 4 $Ast [expr $yb1-$cover] [expr $zb1-$cover] [expr $yb1-$cover] [expr -$zb1+$cover] 
}    
#____________________________________
# BEAM torsional stiffness

set GJbea [expr $Gc*$C250*$beaDepth*pow($beaWidth,3)]
set GAbea [expr $Gc*$beaWidth*$beaDepth*5/6]

uniaxialMaterial Elastic 54 $GJbea
uniaxialMaterial Elastic 55 $GAbea

# Attach torsion to the RC beam section
#                 tag uniTag uniCode secTag
#section Aggregator $secTag $matTag1 $string1 $matTag2 $string2 ....... <-section $sectionTag>
section Aggregator     40    55 Vy      55 Vz    54 T      -section 4 

# #############################################################################################
#COLUMNS#
# Geometry of column elements
#                tag 

geomTransf Linear 1  

# Number of integration points along length of element
set np 5
set nps 10
# Create the columns using Beam-column elements
#               e              tag ndI ndJ nsecs secID transfTag
set eleType forceBeamColumn
element $eleType				1    1   3   $np    10      1
# element $eleType				101  12   3   $np   10      1
element $eleType				2    2   4  $np    10       1
# element $eleType				202  15   4   $np    10     1
# #############################################################################################
#BEAMS#
# Define beam element
# -----------------------------

# Geometry of beam elements
#                tag 
geomTransf Linear 2  
#geomTransf Corotational 3  0 1 0
#geomTransf Corotational 4 -1 0 0
 #1=4
  #2=3
set tranfeloriz 2
set tranfelvert 1
set tranfelorizcor 3
set tranfelvertcor 4 

# Create the NONLINEAR BEAM  using Beam-column elements
#               e              tag ndI ndJ nsecs secID transfTag
set eleType2 forceBeamColumn
# TRAVI
# element $eleType  	    		3   1   2  $np       40   $tranfeloriz
element $eleType2     			4   3   4   $np      40   $tranfeloriz 
#element elasticBeamColumn    4  3    4   8000000  3000000000   1066666666    $tranfeloriz
# Create the ELASTIC BEAM beam element

#                        eleTag iNode jNode  A         E          Iz       transfTag                                          transfTag
 element elasticBeamColumn    11  1    2   8000000  3000000   1066666666    $tranfeloriz

# Create the infill panel element
#                      tag    inode   jnode   -mat  mattag  -dir    dirtag
element zeroLength      5       3       5     -mat   4      -dir     1  

# Define gravity loads
# --------------------

# Set a parameter for the axial load
set P 400.000;                # 10% of axial capacity of columns

# Create a Plain load pattern with a Linear TimeSeries
pattern Plain 1 "Linear" {

    # Create nodal loads at nodes 3 & 4
	#    nd    FX       FY        MZ 
	load  3   0.0   [expr -$P]    0.0
	load  4   0.0   [expr -$P]    0.0
	

}

# ------------------------------
# End of model generation
# ------------------------------

# ------------------------------
# Start of analysis generation
# ------------------------------

# Create the system of equation, a sparse solver with partial pivoting
system BandGeneral

# Create the constraint handler, the transformation method
constraints Transformation

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer RCM

# Create the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
test NormDispIncr 1.0e-12  1000 3

# Create the solution algorithm, a Newton-Raphson algorithm
algorithm Newton

# Create the integration scheme, the LoadControl scheme using steps of 0.1 
integrator LoadControl 0.1

# Create the analysis object
analysis Static

# ------------------------------
# End of analysis generation
# ------------------------------


# ------------------------------
# Perform the gravity analysis
# ------------------------------

# perform the gravity load analysis, requires 10 steps to reach the load level
analyze 10

# Print out the state of nodes 3 and 4
print node 3 4 5

# Print out the state of element 1
print ele 

# ------------------------------
# Perform the cyclic analysis
# ------------------------------

loadConst -time 0.0

set LunitTXT 1
set node 3
set dof  1
set LCol 3125

recorder Node -file node3disp.out -time -node 3 -dof 1 disp;           # disp of node 3 along time
recorder Node -file total_reaction.out  -node 1 2 5  -dof 1 reaction;		# total reaction in direction 1(x)

set IDctrlNode $node;	# node where displacement is read for displacement control
set IDctrlDOF $dof;	# degree of freedom of displacement read for displacement control# characteristics of cyclic analysis	
set iDmax "0.00016 0.00032 0.00064 0.00128 0.00192 0.00256 0.0032 0.00384 0.00512 0.0064 0.00768 0.01024 0.0128 0.01408 0.01664 0.0192 0.02432"

	# vector of displacement-cycle peaks, in terms of storey drift ratio
set Dincr [expr 0.00005*$LCol];	# displacement increment for pushover. you want this to be very small, but not too small to slow down the analysis
		# scale drift ratio by storey height for displacement cycles
set Fact $LCol;						# scale drift ratio by storey height for displacement cycles
set CycleType Full;					# you can do Full / Push / Half cycles with the proc
set Ncycles 1;						# specify the number of cycles at each peak

 # create load pattern for lateral pushover load
 set Hload  1;			# define the lateral load as a proportion of the weight so that the pseudo time equals the lateral-load coefficient when using linear load pattern
  set iPushNode "3 4";			# define nodes where lateral load is applied in static lateral analysis
 pattern Plain 200 Linear {;			# define load pattern -- generalized
	foreach PushNode $iPushNode {
		 load $PushNode $Hload 0.0 0.0
	 }
 }
# ----------- set up analysis parameters
source LibAnalysisStaticParameters.tcl;	# constraintsHandler,DOFnumberer,system-ofequations,convergenceTest,solutionAlgorithm,integrator

#  ---------------------------------    perform Static Cyclic Displacements Analysis
source LibGeneratePeaks.tcl

set fmt1 "%s Cyclic analysis: CtrlNode %.3i, dof %.1i, Disp=%.4f %s";	# format for screen/file output of DONE/PROBLEM analysis
foreach Dmax $iDmax {
	set iDstep [GeneratePeaks $Dmax $Dincr $CycleType $Fact];	# this proc is defined above
	for {set i 1} {$i <= $Ncycles} {incr i 1} {
		set zeroD 0
		set D0 0.0
		foreach Dstep $iDstep {
			set D1 $Dstep
			set Dincr [expr $D1 - $D0]
			integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr
			analysis Static
			# ----------------------------------------------first analyze command------------------------
			set ok [analyze 1]
			# ----------------------------------------------if convergence failure-------------------------
			if {$ok != 0} {
				# if analysis fails, we try some other stuff
				# performance is slower inside this loop	global maxNumIterStatic;	    # max no. of iterations performed before "failure to converge" is ret'd
				# if {$ok != 0} {
					# puts "Trying Newton with Initial Tangent .."
					# test NormDispIncr   1.0e-4 2000 0
					# algorithm Newton -initial
					# set ok [analyze 1]
					# test $testTypeStatic $TolStatic      $maxNumIterStatic    0
					# algorithm $algorithmTypeStatic
				# }
				# if {$ok != 0} {
					# puts "Trying Broyden .."
					# algorithm Broyden 8
					# set ok [analyze 1 ]
					# algorithm $algorithmTypeStatic
				# }
				if {$ok != 0} {
					puts "Trying NewtonWithLineSearch .."
					algorithm NewtonLineSearch 0.8 
					set ok [analyze 1]
					algorithm $algorithmTypeStatic
				}
				if {$ok != 0} {
					set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
					puts $putout
					return -1
				}; # end if
			}; # end if
			# -----------------------------------------------------------------------------------------------------
			set D0 $D1;			# move to next step
		}; # end Dstep
	};		# end i
};	# end of iDmaxCycl
# -----------------------------------------------------------------------------------------------------
if {$ok != 0 } {
	puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
} else {
	puts [format $fmt1 "DONE"  $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
}
