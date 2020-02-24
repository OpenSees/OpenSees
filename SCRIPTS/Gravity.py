

# ------------------------------
# Gravity load analysis
# ------------------------------

# Generate the model and load input variables
ops.source('MVLEM_3d_Model.tcl')
ops.source('SFI_MVLEM_3d_Model.tcl')

# Create a Plain load pattern with a linear TimeSeries
ops.pattern('Plain',1,'"Linear"')
 
    # Create the nodal load - command: load nodeID xForce yForce
ops.load(81,0.0,'[expr','-$N/2.0]',0.0,0.0,0.0,0.0)
ops.load('$IDctrlNode',0.0,'[expr','-$N/2.0]',0.0,0.0,0.0,0.0)

# ------------------------------
# Analysis generation
# ------------------------------

# Create the integration scheme, the LoadControl scheme using steps of 0.1 
ops.integrator('LoadControl',0.01)
 
# Create the system of equation, a sparse solver with partial pivoting
ops.system('BandGeneral')

# Create the convergence test, the norm of the residual with a tolerance of 
# 1e-5 and a max number of iterations of 100 
ops.test('NormDispIncr','$Tol',100,0)

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
ops.numberer('RCM')

# Create the constraint handler, the transformation method
ops.constraints('Transformation')

# Create the solution algorithm, a Newton-Raphson algorithm
ops.algorithm('Newton','#','-initialThenCurrent')

# Create the analysis object
ops.analysis('Static')

# Run analysis
ops.analyze(100)
