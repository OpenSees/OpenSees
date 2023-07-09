#================================================================
# CreateMaterialModels.tcl
# --Create the all elements from file elements.dat
# Zhaohui Yang
# April 6, 2002
#================================================================

#------------------------------------------------------
# Drucker-Prager model for sand
#------------------------------------------------------
#Yield surface 
set DPys "-DP"

#potential surface
set DPps "-DP 0.0"

# Scalar evolution law: linear hardening coef = 1.0
set ES1  "-Leq  0.0"

# Tensorial evolution law: linear hardening coef. = 0.0
set ET1  "-Linear  0.0"

# initial stress
set stressp "-0.01 0 0  0 -0.01 0  0 0 -0.01"

# EPState
set E  [expr 17400*$kPa]; #Young's Modulus or any number
set Eo [expr 17400*$kPa]; #Young's Modulus at reference pressure (1 P_atm)
set Vs 0.35;  #Poisson's Ratio

set rho_sand [expr 1478.0*$kg/pow($m, 3)]; #unit weight
set phi [expr 37*$deg ]; #Friction angle
set alpha  [expr 2.0*sin($phi)/( 1.7321*(3-sin($phi)) )]

#_____ ___E___Eo__v______rho__________________alpha___k
set EPS "$E  $Eo $Vs $rho_sand -NOD 1 -NOS 2  $alpha  0.0  -stressp $stressp"
#
# where
#alpha = 2 sin(phi) / (3^0.5) / (3-sin(phi) ), phi is the friction angle
# and k is the cohesion

# Create nDMaterial using Template Elastic-PLastic Model
#nDMaterial Template3Dep 1 -YS $DPys -PS $DPps -EPS $EPS -ELS1 $ES1 -ELT1 $ET1
nDMaterial Template3Dep 1 -YS $DPys -PS $DPps -EPS $EPS -ELS1 $ES1 -ELT1 $ET1
