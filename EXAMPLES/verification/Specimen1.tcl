#______________________________________________________________________ 
# Example Analysis of
# Sezen Specimen Column-1

# CREATION OF MODEL
		model BasicBuilder -ndm 2 -ndf 3
		puts "Cyclic Analysis Sezen Column"
		set dataDir Results
		file mkdir $dataDir

# CREATION OF NODES
#    		  tag	X   Y
		node  100  0.0 0.0
		node  200  0.0 116.0

#  			 tag   DX DY RZ
		fix  100   1  1  1
		fix  200   0  0  1

# NODE FOR SPRING
		node  101  0.0 0.0
		node  102  0.0 116.0


# Column Section
#  MATERIALS AND SECTION
set Conctag		    1
set Covertag		2
set Steeltag		3 
set Ec [expr sqrt(3350)*57000]



# Concrete Materials for confined and unconfined concrete

      uniaxialMaterial Concrete01  $Covertag    -3.06    -0.002  -3.06  -0.007
#     uniaxialMaterial Concrete04  $Covertag    -3.06    -0.002  -0.006 [expr $Ec/1000 ]
      uniaxialMaterial Concrete04  $Conctag     -3.35    -0.003  -0.0185 [expr $Ec/1000 ]


# Materials for Reinforrcing Steel.

#	uniaxialMaterial ReinforcingSteel $Steeltag 59 93.5 29000 700 0.002 0.15
	uniaxialMaterial Steel01          $Steeltag 59 29000.0  0.02

set ColSec 		1 
set colWidth 	18.0
set colDepth 	18.0
set cover  		2.0; 	#concrete cover
set As    		1.0; 		#area of longitudinal bars

# some variables derived from the parameters
set y1 [expr $colDepth/2.0]
set z1 [expr $colWidth/2.0]

section Fiber $ColSec  {

    # Create the concrete core fibers
    patch rect $Conctag 10 1 [expr $cover-$y1] [expr $cover-$z1] [expr $y1-$cover] [expr $z1-$cover]

    # Create the concrete cover fibers (top, bottom, left, right)
    patch rect $Covertag 10 1  [expr -$y1] [expr $z1-$cover] $y1 $z1
    patch rect $Covertag 10 1  [expr -$y1] [expr -$z1] $y1 [expr $cover-$z1]
    patch rect $Covertag 2 1  [expr -$y1] [expr $cover-$z1] [expr $cover-$y1] [expr $z1-$cover]
    patch rect $Covertag 2 1  [expr $y1-$cover] [expr $cover-$z1] $y1 [expr $z1-$cover]
	
	    # Define the reinforcement explicitly using fiber command
    #     yloc  zloc  area matTag
    fiber  6.5625 -6.5625 1   $Steeltag  0 
    fiber  6.5625  0      1   $Steeltag  0 
    fiber  6.5625  6.5625 1   $Steeltag  0 
    fiber  0.0     6.5625 1   $Steeltag  0 
    fiber  0.0    -6.5625 1   $Steeltag  0
    fiber -6.5625 -6.5625 1   $Steeltag  0
    fiber -6.5625  0.0    1   $Steeltag  0
    fiber -6.5625  6.5625 1   $Steeltag  0
	
	
	

} 

 set Vcr  23.418982
 set Vpeak  69.403443
 set d_cr  0.006915
 set d_p  0.114425
 set d_u  0.4
 set d_alf  6.144612
 set V1  2380.008376
 set D1  0.000785
 set V2  3881.164156
 set D2  0.002882
 set V3  4376.958758
 set D3  0.003452

 
# Springs 
set ShearTag 		4
set rigidMatTag 	5
set SlipTag		    6
set Axial				7
#  MATERIALS FOR SPRING
#________MATERIALS____________________________________________________________________________________________
	set E		[expr ($Vcr)/($d_cr)]
	set ecr		$d_cr
	set E2		[expr ($Vpeak - $Vcr)/($d_p - $d_cr)]
	set eyp		$d_p
	set eye		$d_u
	set ealf	$d_alf
#_____________________________________________________________________________________________________________
uniaxialMaterial Shear01 $ShearTag $E $ecr $E2 $eyp $eye $ealf $E [expr -$ecr] $E2 [expr -$eyp] [expr -$eye] [expr -$ealf]   
uniaxialMaterial Elastic $rigidMatTag 9.9e10
uniaxialMaterial Axial01 $Axial 101 100 1.08e4 $d_u $d_alf

# this part is to use hyteretic material for Slip.
		set px 1
		set py 1
		set dam1 0.0
		set dam2 0.0
		set beta 0.5
#

		uniaxialMaterial Hysteretic $SlipTag $V1 $D1 $V2 $D2 $V3 $D3\
					[expr -$V1] [expr -$D1] [expr -$V2] [expr -$D2] [expr -$V3] [expr -$D3] \
					$px $py $dam1 $dam2 $beta
#_____________________________________________________________________________________________________________


# FLEXURAL SECTION
geomTransf Linear 1
set ite 5
set element forceBeamColumn
element $element 1 101 102 $ite $ColSec 1 

#element forceBeamColumn 1 101 102 1 "HingeRadau $ColSec 27 $ColSec 27 $ColSec" 
#		Axial, Slip and Shear  Spring for top and bottom

	element zeroLength 2 100 101 -mat $ShearTag $rigidMatTag $SlipTag -dir 1 2 3
	element zeroLength 3 102 200 -mat $rigidMatTag $Axial $SlipTag -dir 1 2 3 






#________________________________________________________________________________________
#  RECORDER
recorder Node    -file $dataDir/RHS1.out 	-node 100 -dof 1  reaction
recorder Node    -file $dataDir/DHS1.out 	-node 200 -dof 1  disp
recorder Node    -file $dataDir/AHS1.out       	-node 200 -dof 2  disp


#_________________________________________________________________________________________
# Initial axial load 
set P -150.0 ;#Compression axial load
# Constant gravity load
pattern Plain 1 Constant {load 200 0.0 $P 0.0 }
#___Analysis options for gravity loads
		system ProfileSPD
		constraints Penalty 1.0e12 1.0e12
		numberer RCM
		test NormDispIncr 1.0e-8 25 1
		algorithm KrylovNewton
		integrator LoadControl 0 1 0 0
		analysis Static
		analyze 1
		wipeAnalysis

# 
set nSteps 9457; #Final point of analysis how many points....
set dlamda 0.0001; #the inverval
set er 1.0e-3;	

timeSeries Path 1 -dt $dlamda -filePath Specimen1.acc;
pattern Plain 2 1 { sp 200 1 1 }

###############################################
# Set Analysis options for transverse loading
###############################################
# Create the system of equation
system BandGeneral
# Create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer RCM
# Create the integration scheme, 
#LoadControl scheme using constant steps of dlamda
#                      dlamda1 Jd minLamda maxLamda
		integrator LoadControl $dlamda 1  $dlamda $dlamda
#      integrator DisplacementControl 200     1  $dlamda 1  $dlamda   $dlamda

constraints Penalty 1.0e14  1.0e14

#########################
# Create the convergence test, Norm of the displacement increment
#                 tol     maxNumIter printFlag
test NormDispIncr $er 100         5
# Create the solution algorithm, a Newton-Rahpson algorithm is created
#algorithm Newton
algorithm KrylovNewton
# create the analysis object 
analysis Static 
#########################################
# Analyze model with transverse loading #
#########################################
set ok 0
set n 1
while {$n < $nSteps && $ok == 0} {
	set ok [analyze 1]
	if {$ok != 0} {
		test NormDispIncr $er 1000 5     ;# increase max number of iterations 
		algorithm ModifiedNewton -initial      ;# use initial elastic stiffness for NewtonRaphson
		#puts "Time Newton Initial [getTime]" ;# output time algorithm was changed
		set ok [analyze 1 ] 			;# analyse 1 step with new settings
		#set ans [gets stdin] 		        ;# pauses tcl script
		algorithm Newton                        ;# restore algorithm to Newton with current stiffness
		test NormDispIncr $er 100 1          ;# restore max number of iterations
	}
	incr n
puts $n
}

if {$ok != 0} {
	puts "ANALYSIS FAILED"
} else {
	puts "SUCCESSFULL"
}

wipe




puts "end of analysis"




