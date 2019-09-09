# ************************************************************************************************************
#              ***************************************************************************
#                                  *************************************
#				Mechanical Characterization of Integrally-Attached Timber Plate Structures
#                                   Aryan Rezaie Rad, EPFL, Suisse, 2018                  
#                                  *************************************
#              ***************************************************************************
# ************************************************************************************************************
#                           Section Name : In-Plane, Rectangular, Isotropic

# In-Plane Action of Forces: Timoshenko Formulation

# ------------------------------ Load Control Analysis -------------------------------------------------------------------

set Dmax 0.2
set Dincr [expr 0.005*$Dmax ]
set Nsteps [expr int($Dmax/$Dincr)];

constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-4 1000
algorithm Newton
analysis Static

integrator DisplacementControl  23 1 $Dincr
set ok [analyze $Nsteps] 

if {$ok != 0 } {
	puts " analysis diverged at displacement [nodeDisp 23 1] "
} else {
	puts " DONE   [nodeDisp 23 1] "
}
