#########################################################################
#---------------Department of Civil, Structural and Environmental Engineering--------------------------#
#----------------------------------------University at Buffalo-----------------------------------------#
#Modeling of elastomeric bearings in compression and tension for response history analysis#
#Written By: Manish Kumar (mkumar2@buffalo.edu)#
#Date: July, 2011 #
#########################################################################

#Units: kips, in, sec

#Remove existing model
wipe
wipeAnalysis

#-----------------------------------------------------------
#User Defined Parameters
#-----------------------------------------------------------
set g 9.810;                 	# Acceleration due to gravity
set pi 3.14159;            		# Value of pi
set p_axial 4.3e+06;            # Axial pressure
set G 0.870e+06;				# Effective shear modulus
set K 2000e+06;                 # Bulk modulus of rubber
set ts 1.016e-03;               # Thickness of steel shim plates
set tr 2.286e-03;             	# Thickness of a single rubber layer
set n 20;                      	# Number of rubber layers
set D1 0.1e-03;               	# Internal diameter of lead rubber bearing
set D2 176e-03;               	# Outer diameter of lead rubber bearing
set tc 7.874e-03;               # Bearing cover
set	a1	234856.6921;
set	a2	-4777446.405;
set	a3	252445116.1;
set	b1	3327.2536;
set	b2	1399628.929;
set	b3	74.29133858;
set	c1	103.740365;
set	c2	360.0400902;
set	c3	0.8743;
set	c4	1159.451138;
set kc 20.0;					# Cavitation parameter
set PhiM 0.5;					# Maximum damage index
set ac 1.0;						# Strength reduction parameter
set DampingRatio 3;				# Damping ratio in percentage for calculation of equivalent yield strength in horizontal
set Strain 100;					# Strain used for calculation of equivalent yield strength in horizontal

#-------------------------------------------------------------------
#Derived Parameters
#-------------------------------------------------------------------

set A [expr ($pi/4)*(($D2+$tc)*($D2+$tc)-$D1*$D1)];             	# Bonded area
set AL [expr $pi*$D1*$D1/4];										# Internal hole/lead area
set S [expr ($D2-$D1)/(4*$tr)];                  					# Shape factor
set Tr [expr $n*$tr];                       						# Total rubber thickness
set h [expr $Tr + ($n-1)*$ts];          							# Total height of bearing
set r [expr $D2/$D1];                   							# Outer to inner diameter ratio
set M [expr $p_axial*$A/$g];										# Equivalent mass for period calculation
if {$D1 == 0} {														# Diameter modification factor
   set F 1.0
} else {
set F [expr ($r*$r+1)/(($r-1)*($r-1)) + (1+$r)/((1-$r)*log($r))];   
}

#For horizontal motion
set uy_h 0.007;                   									# Yield displacement of elastomeric bearing in horizontal direction
set LateralDisp [expr 0.01*$Strain*$Tr];							# Displacement used for calculation of yield strength
set Fy_h [expr 0.5*$pi*0.01*$DampingRatio*$G*$A*$LateralDisp/$Tr];  # Yield strength of elastomeric bearing in horizontal direction
set k1 [expr $Fy_h/$uy_h];             								# Elastic stiffness of bearing
set k2 [expr $G*$A/$Tr];                							# Post-yield stiffness of bearing
set Ccr [expr 2*sqrt(($k2)*$M)];									# Critical damping
set cd [expr 0.01*$DampingRatio*$Ccr];  							# Damping in the elastomeric bearing

#For vertical motion
set Ec [expr 1.0/((1.0/(6.0*$G*$S*$S*$F))+(4.0/(3.0*$K)))];      	# Compressive modulus of elastomeric bearing
set E [expr 3*$G];                                                  # Elastic modulus
set I [expr ($pi/64)*(($D2+$tc)**4-$D1**4)];                        # Moment of inertia of bearing
set rg [expr sqrt($I/$A)];                                          # Radius of gyration 
set Kpre [expr $A*$Ec/$Tr ];                                       	# Pre-cavitation stiffness in Tension
set Kpost [expr $A*$E/$Tr ];                                       	# Post-cavitation stiffness in Tension
set Fc [expr 3*$G*$A];                                              # Cavitation force
set uc [expr $Fc/$Kpre ];                                           # Cavitation displacement
set Er [expr $Ec/3];                                     			# Rotation modulus of bearing
set As [expr $A*$h/$Tr];                                 			# Adjusted shear area of bearing
set Is [expr $I*$h/$Tr];                                			# Adjusted moment of inertia of bearing
set Pe [expr $pi*$pi*$Er*$Is/($h*$h)];                   			# Euler buckling load of bearing
set Pcr [expr -sqrt($Pe*$G*$As)];     								# Critical buckling load in compression
set ucr [expr $Pcr/$Kpre];                                          # Critical displacement in compression


#print parameters
#puts " A : $A, Tr: $Tr, h: $h, S: $S, F: $F"
#puts " Horizontal motion: k1: $k1, k2: $k2, Fy_h: $Fy_h, uy_h: $uy_h, cd: $cd"
#puts " Vertical motion: Ec: $Ec, E: $E, r: $r, Kpre: $Kpre, Kpost: $Kpost, Fc: $Fc, uc: $uc, Pcr: $Pcr, ucr: $ucr"

#---------------------------------------------------------------------------------------------
#Start of model generation
#---------------------------------------------------------------------------------------------

# Elastomeric bearing is modeled as 2 node and 3 DOF element of height h
#Create Model Builder
model basic -ndm 3 -ndf 6

# Create nodes
node 1 0 0 0
node 2 0 $h 0

# Define single point constraints (Constrain the isolator against rotation at both nodes)
fix 1 1 1 1 1 1 1
fix 2 0 0 0 1 1 1

#Define material and element for elastomeric bearing
element HDR 1 1 2 $Fy_h $uy_h $G $K $D1 $D2 $ts $tr $n $a1 $a2 $a3 $b1 $b2 $b3 $c1 $c2 $c3 $c4 0 1 0 1 0 0 10 0.75 1.0 0.5 0.0 0
#-------------------------------------------------------------------------------------------------
#Define loads
#-------------------------------------------------------------------------------------------------

# Apply gravity load on isolated mass
set P [expr $p_axial*$A]

#Create a plain load pattern with linear timeseries

pattern Plain 1 "Linear" {
	
	load 2 0.0 [expr -$P] 0.0 0.0 0.0 0.0
	sp 2 1 [expr $Tr]
}

#---------------------------------------------------------------------------
#Start of analysis generation
#---------------------------------------------------------------------------
system BandSPD
constraints Transformation
numberer RCM
test NormDispIncr 1.0e-15 10 3
algorithm Newton
integrator LoadControl 1
analysis Static

# -------------------------------------------------------------------------
# Start of recorder generation
# -------------------------------------------------------------------------

# Create a recorder to monitor nodal displacements and force in elastomeric bearing element
#recorder Node -file dispGravity.out -time -node 2 -dof 1 2 3 disp
#recorder Element -file forceGravity.out -time -ele 1 force

# -------------------------------------------------------------------------
# Finally perform the analysis
# -------------------------------------------------------------------------

set ok 0

set ok [analyze 1]

#Print a message to indicate if analysis was successful or not
puts [nodeDisp 2 1]
puts [nodeDisp 2 2]
puts $Tr
if {([nodeDisp 2 1] == $Tr) && ([nodeDisp 2 2] == -0.00027115536497439714)} { puts "SUCCESS"} else {puts "FAILURE"}
