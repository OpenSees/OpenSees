from openseespy.opensees import *

# Remove existing model
wipe()

# EXAMPLE 3 (Kim and Constantinou, 2023 https://doi.org/10.1002/eqe.3797)
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
model('basic', '-ndm', 3, '-ndf', 6)

# Create nodes
node(1, 0, 0, 0)  # End i
node(2, 0, 0, 0)  # End j

# Define single point constraints
fix(1, 1, 1, 1, 1, 1, 1)

# Define friction models
tagTemp = 1
tagVel = 0
tagPres = 0
velRate = 100
kTmodel = 1  # kT = 1/2 at 200 degree celsius

#----------------------------------------------------------------------------
# Bring material models and define element
#----------------------------------------------------------------------------

# Creating material for compression and rotation behaviors
uniaxialMaterial('Elastic', 1, kvc)
uniaxialMaterial('Elastic', 2, 10.)

tagT = 1

# Define TripleFrictionPendulumX element
element('TripleFrictionPendulumX', 1, 1, 2, tagT, 1, 2, 2, 2, tagPres, tagTemp, tagVel, mu1, mu2, mu3,
            L1, L2, L3, d1, d2, d3, b1, b2, b3, P, uy, kvt, minFv, tol, Pref1, Pref2, Pref3, Diffu, Conduct,
            Temperature0, velRate, kTmodel, 1)

#----------------------------------------------------------------------------
# Apply gravity load
#----------------------------------------------------------------------------

# Create a plain load pattern with linear timeseries
timeSeries('Linear', 1)
pattern('Plain', 1, 1)
load(2, 0.0, 0.0, -P, 0.0, 0.0, 0.0)

#----------------------------------------------------------------------------
# Start of analysis generation (Gravity)
#----------------------------------------------------------------------------

system('BandSPD')
constraints('Transformation')
numberer('RCM')
test('NormDispIncr', 1.0e-15, 10)
algorithm('Newton')
integrator('LoadControl', 1)
analysis('Static')

#----------------------------------------------------------------------------
# Analysis (Gravity)
#----------------------------------------------------------------------------

analyze(1)
print("Gravity analysis completed SUCCESSFULLY")

#----------------------------------------------------------------------------
# Start of analysis generation
# (Sinusoidal; Two cycles of 5s period and 508mm amplitude)
#----------------------------------------------------------------------------

loadConst('-time', 0.0)

# analysis time step
dt = 0.008

# excitation time step
dt1 = 0.001

timeSeries('Trig', 11, dt, 10, 5, '-factor', 0.508, '-shift', 0, '-zeroShift', 0)
pattern('MultipleSupport', 2)
groundMotion(1,'Plain', '-disp', 11)
imposedMotion(2,2,1)

#----------------------------------------------------------------------------
# Start of recorder generation (Sinusoidal)
#----------------------------------------------------------------------------

# Set up recorder
OutDir = "EXAMPLE3"  # Output folder
OutFile1 = "TEMPERATURE.out"
OutFile2 = "DISP.out"
OutFile3 = "FORCE.out"
OutFile4 = "COMPDISP.out"

import os
os.makedirs(OutDir, exist_ok=True)

recorder('Element', '-file', os.path.join(OutDir, OutFile1), '-time', '-ele', 1, 'Parameters')
recorder('Node', '-file', os.path.join(OutDir, OutFile2), '-time', '-nodes', 2, '-dof', 1, 2, 3, 'disp')
recorder('Element', '-file', os.path.join(OutDir, OutFile3), '-time', '-ele', 1, 'basicForce')
recorder('Element', '-file', os.path.join(OutDir, OutFile4), '-time', '-ele', 1, 'compDisplacement')

#----------------------------------------------------------------------------
# Analysis (Sinusoidal)
#----------------------------------------------------------------------------

wipeAnalysis()
system('SuperLU')
constraints('Transformation')
numberer('Plain')
integrator('Newmark', 0.5, 0.25)
test('NormDispIncr', 1.0e-5, 20)
algorithm('Newton')
analysis('Transient')

# set some variables
tFinal = 10
tCurrent = getTime()
ok = 0

# Perform the transient analysis
while ok == 0 and tCurrent < tFinal:
    ok = analyze(1, dt)

    # if the analysis fails, try initial tangent iteration
    if ok != 0:
        print("Regular Newton failed... Let's try an initial stiffness for this step")
        test('NormDispIncr', 1.0e-12, 100)
        algorithm('ModifiedNewton', '-initial')
        ok = analyze(1, dt)
        if ok == 0:
            print("That worked... Back to regular Newton")
        test('NormDispIncr', 1.0e-12, 10)
        algorithm('Newton')

    tCurrent = getTime()

# Print a message to indicate if analysis was successful or not
if ok == 0:
    print("Transient analysis completed SUCCESSFULLY")
else:
    print("Transient analysis completed FAILED")
