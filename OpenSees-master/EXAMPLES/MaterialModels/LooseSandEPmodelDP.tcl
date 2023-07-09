# 
# Template elasto-plastic material model for loose sand, Drucker_Prager yield surface
# rho = 1.7ton/m^3 
# Young's modulus = 148500 kpa, 
# Poisson's ratio = 0.335
# Friction Angle = 29deg, 
# Cohesion k = 0.0

# Aug. 12, 2002  Feng Xiong
# 														 Boris Jeremic (@ucdavis.edu)


# Yield surface
set DPys "-DP"

# Potential surface
set DPps "-DP 0.1"

# Scalar evolution law: linear hardening coef=1.0
set ES1 "-Leq 1.10"

# Initial stress
set stressp 0.1 0 0 0 0.1 0 0 0 0.1

# EPState
#------------E--------Eo--------v----rho-----------------alpha----k
set EPS " 148500.0  148500.0  0.335  1.7  -NOD 0 -Nos 2   0.223  0.0  -stressp $stressp"

# Create nDMaterial using Template Elastic-Plastic Model
nDMaterial Template3Dep 550 -YS $DPys -PS $DPps -EPS $EPS -ELS1 $ES1

 
