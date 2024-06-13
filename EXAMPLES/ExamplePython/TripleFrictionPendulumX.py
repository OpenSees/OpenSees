import sys
#sys.path.append(r'C:\OpenSees-3.3.0\Win64\bin')
#import opensees as op

import openseespy.opensees as op 

# Units: N, m, sec
# Remove existing model
op.wipe()

# TFP Geomoetry of Configuration A (Kim and Constantinou, 2023 https://doi.org/10.1002/eqe.3797)
#----------------------------------------------------------------------------
# User Defined Parameters
#----------------------------------------------------------------------------

# TFP Geomoetry of Configuration A
L1 = 0.3937  # Effective radii (m)
L2 = 3.7465
L3 = 3.7465
d1 = 0.0716  # Actual displacement capacity (m)
d2 = 0.5043
d3 = 0.5043
b1 = 0.508  # Diameter of the rigid slider and the two inner slide plate (m)
b2 = 0.711
b3 = 0.711
r1 = b1 / 2  # Radius of the rigid slider and the two inner slide plate (m)
r2 = b2 / 2
r3 = b3 / 2
Thickness2 = 0.02  # Thickness of concave plate (m)
Thickness3 = 0.02


uy = 0.001  # Yield displacement (m)
kvc = 8000000000.  # Vertical compression stiffness (N/m)
kvt = 1.  # Vertical tension stiffness (N/m)
minFv = 0.1  # Minimum compression force in the bearing (N)

g = 9.81  # Gravity acceleration (m/s^2)
P = 13345e+03  # Load on top of TFP
Mass = P / g  # Mass on top of TFP
tol = 1.e-5  # Relative tolerance for checking convergence

# Heat parameters
Diffu = 0.444e-5  # Thermal diffusivity (m^2/sec)
Conduct = 18  # Thermal conductivity (W/m*Celsius)
Temperature0 = 20  # Initial temperature (Celsius)
tagT2 = 2 					  # 1 = indefinite plate thickness / 2 = finite plate thickness

# Friction coefficients (reference)
mu1 = 0.01
mu2 = 0.04
mu3 = 0.08

# Reference Pressure
Pref1 = P / (r1 * r1 * 3.141592)
Pref2 = P / (r2 * r2 * 3.141592)
Pref3 = P / (r3 * r3 * 3.141592)

#----------------------------------------------------------------------------
# Start of model generation
#----------------------------------------------------------------------------

# Create Model Builder
op.model('basic', '-ndm', 3, '-ndf', 6)

# Create nodes
op.node(1, 0, 0, 0)  # End i
op.node(2, 0, 0, 0)  # End j

# Define single point constraints
op.fix(1, 1, 1, 1, 1, 1, 1)

# Define friction models
tagTemp = 1
tagVel = 1
tagPres = 0
velRate = 100
kTmodel = 1  # kT = 1/2 at 200 degree celsius

#----------------------------------------------------------------------------
# Bring material models and define element
#----------------------------------------------------------------------------

# Creating material for compression and rotation behaviors
op.uniaxialMaterial('Elastic', 1, kvc)
op.uniaxialMaterial('Elastic', 2, 10.)

tagT = 1

# Define TripleFrictionPendulumX element
op.element('TripleFrictionPendulumX', 1, 1, 2, tagT, tagT2, 1, 2, 2, 2, tagPres, tagTemp, tagVel, mu1, mu2, mu3,
            L1, L2, L3, d1, d2, d3, b1, b2, b3, Thickness2, Thickness3, P, uy, kvt, minFv, tol, Pref1, Pref2, Pref3, Diffu, Conduct,
            Temperature0, velRate, kTmodel, 1)

#----------------------------------------------------------------------------
# Apply gravity load
#----------------------------------------------------------------------------

# Create a plain load pattern with linear timeseries
op.timeSeries('Linear', 1)
op.pattern('Plain', 1, 1)
op.load(2, 0.0, 0.0, -P, 0.0, 0.0, 0.0)

#----------------------------------------------------------------------------
# Start of analysis generation (Gravity)
#----------------------------------------------------------------------------

op.system('BandSPD')
op.constraints('Transformation')
op.numberer('RCM')
op.test('NormDispIncr', 1.0e-15, 10)
op.algorithm('Newton')
op.integrator('LoadControl', 0.1)
op.analysis('Static')

#----------------------------------------------------------------------------
# Analysis (Gravity)
#----------------------------------------------------------------------------

op.analyze(10)
print("Gravity analysis completed SUCCESSFULLY")

#----------------------------------------------------------------------------
# Start of analysis generation
# (Sinusoidal; Ten cycles of 5s period and 600mm amplitude)
#----------------------------------------------------------------------------

op.loadConst('-time', 0.0)

# analysis time step
dt = 0.008

# excitation time step
dt1 = 0.001

op.timeSeries('Trig', 11, dt, 50, 5, '-factor', 0.6, '-shift', 0, '-zeroShift', 0)
op.pattern('MultipleSupport', 2)
op.groundMotion(1,'Plain', '-disp', 11)
op.imposedMotion(2,2,1)

#----------------------------------------------------------------------------
# Start of recorder generation (Sinusoidal)
#----------------------------------------------------------------------------

# Set up recorder
OutDir = "EXAMPLE_PY"  # Output folder
OutFile1 = "TEMPERATURE_FINITE_DEPTH.out"
OutFile2 = "DISP_FINITE_DEPTH.out"
OutFile3 = "FORCE_FINITE_DEPTH.out"
OutFile4 = "COMPDISP_FINITE_DEPTH.out"

# OutFile1 = "TEMPERATURE_INDEFINITE_DEPTH.out"
# OutFile2 = "DISP_INDEFINITE_DEPTH.out"
# OutFile3 = "FORCE_INDEFINITE_DEPTH.out"
# OutFile4 = "COMPDISP_INDEFINITE_DEPTH.out"

import os
os.makedirs(OutDir, exist_ok=True)

op.recorder('Element', '-file', os.path.join(OutDir, OutFile1), '-time', '-ele', 1, 'Parameters')
op.recorder('Node', '-file', os.path.join(OutDir, OutFile2), '-time', '-nodes', 2, '-dof', 1, 2, 3, 'disp')
op.recorder('Element', '-file', os.path.join(OutDir, OutFile3), '-time', '-ele', 1, 'basicForce')
op.recorder('Element', '-file', os.path.join(OutDir, OutFile4), '-time', '-ele', 1, 'compDisplacement')

#----------------------------------------------------------------------------
# Analysis (Sinusoidal)
#----------------------------------------------------------------------------

op.wipeAnalysis()
op.system('SuperLU')
op.constraints('Transformation')
op.numberer('Plain')
op.integrator('Newmark', 0.5, 0.25)
op.test('NormDispIncr', 1.0e-5, 20)
op.algorithm('Newton')
op.analysis('Transient')

# set some variables
tFinal = 50
tCurrent = op.getTime()
ok = 0

# Perform the transient analysis
while ok == 0 and tCurrent < tFinal:
    ok = op.analyze(1, dt)
    tCurrent = op.getTime()

# Print a message to indicate if analysis was successful or not
if ok == 0:
    print("Transient analysis completed SUCCESSFULLY")
else:
    print("Transient analysis completed FAILED")

op.wipe()
