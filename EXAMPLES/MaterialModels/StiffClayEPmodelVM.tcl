# 
# Template elasto-plastic material model for stiff clay, von Mises yield surface
# rho = 1.8 ton/m^3  
# Young's modulu = 423000 kpa, 
# Possion Ratio = 0.406  
# cohesion k (=f(Cu)) = 75kpa  

# Aug. 12, 2002  Feng Xiong


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
set EPS "423000.0  423000.0   0.406  1.8  -NOD 0 -Nos 1  75 -stressp $stressp"

# Create nDMaterial using Template Elastic-Plastic Model
nDMaterial Template3Dep 1200 -YS $DPys -PS $DPps -EPS $EPS -ELS1 $ES1

 
