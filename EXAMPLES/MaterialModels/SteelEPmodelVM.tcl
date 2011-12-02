# 
# Template elasto-plastic material model for steel, von Mises yield surface
# rho = 1.7 ton/m^3
# Young's modulus = 210000 kpa, 
# Possion ratio = 0.28 
# shear strenght k (=f(Cu)) = 20kpa  

# Aug. 12, 2002  Feng Xiong
# 														 Boris Jeremic (@ucdavis.edu)


# Yield surface
set DPys "-VM"

# Potential surface
set DPps "-VM"

# Scalar evolution law: linear hardening coef=1.0
set ES1 "-Leq 1.10"

# Initial stress
set stressp 0.1 0 0 0 0.1 0 0 0 0.1

# EPState
#------------E--------Eo--------v----rho---------------k=f(Cu)
set EPS "210000.0  210000.0  0.28  1.7  -NOD 0 -Nos 1 20 -stressp $stressp"

# Create nDMaterial using Template Elastic-Plastic Model
nDMaterial Template3Dep 1600 -YS $DPys -PS $DPps -EPS $EPS -ELS1 $ES1

 
