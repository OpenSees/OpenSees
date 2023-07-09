# 
# Template elasto-plastic material model for medium dense sand, Drucker_Prager yield surface
# rho = 2.0 ton/m^3 
# Young's modulus = 270000 kpa, 
# Possion Ratio = 0.35
# Friction Angle = 37deg,  
# kappa = 0.06

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
set EPS " 270000.0  270000.0  0.350  2.0  -NOD 0 -Nos 2   0.290  0.0  -stressp $stressp"

# Create nDMaterial using Template Elastic-Plastic Model
nDMaterial Template3Dep 750 -YS $DPys -PS $DPps -EPS $EPS -ELS1 $ES1

 
