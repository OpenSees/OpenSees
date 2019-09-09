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

constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-6 1000
algorithm Newton
integrator LoadControl 0.05
analysis Static
analyze 20
