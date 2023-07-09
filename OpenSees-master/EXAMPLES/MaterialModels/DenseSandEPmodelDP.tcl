# 
# Template elasto-plastic material model for dense sand, Drucker_Prager yield surface
# Friction Angle = 40deg, 
# rho = 2.1ton/m^3  
# Young's modulus = 351000 kpa, 
# Poisson's ratio = 0.35

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
set EPS " 351000.0  351000.0  0.350  2.1  -NOD 0 -Nos 2   0.315  0.0  -stressp $stressp"

# Create nDMaterial using Template Elastic-Plastic Model
nDMaterial Template3Dep 850 -YS $DPys -PS $DPps -EPS $EPS -ELS1 $ES1

 
