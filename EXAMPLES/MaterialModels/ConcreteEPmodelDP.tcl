# 
# Template elasto-plastic material model for concrete, Drucker-Prager yield surface
# Young's modulus E=30000 kpa, 
# Poisson's ratio v=0.18, 
# Friction angle = 40deg
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
set EPS " 30000.0  30000.0     0.18  2.4  -NOD 0 -Nos 2   0.315  0.0  -stressp $stressp"

# Create nDMaterial using Template Elastic-Plastic Model
nDMaterial Template3Dep 1500 -YS $DPys -PS $DPps -EPS $EPS -ELS1 $ES1

 
