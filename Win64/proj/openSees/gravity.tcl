# ------------------------------
# Gravity load analysis
# ------------------------------

# Generate the model and load input variables
source SFI_MVLEM_SP4.tcl

# Create a Plain load pattern with a linear TimeSeries
pattern Plain 1 "Linear" {
 
    # Create the nodal load - command: load nodeID xForce yForce
    load $IDctrlNode  0.0 [expr  -$N]  0.0     
}			

# ------------------------------
# Analysis generation
# ------------------------------

# Create the integration scheme, the LoadControl scheme using steps of 0.1 
integrator LoadControl 0.1
 
# Create the system of equation, a sparse solver with partial pivoting
system BandGeneral		

# Create the convergence test, the norm of the residual with a tolerance of 
# 1e-5 and a max number of iterations of 100 
test NormDispIncr $Tol  100 0

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer RCM	

# Create the constraint handler, the transformation method
constraints Transformation	

# Create the solution algorithm, a Newton-Raphson algorithm
algorithm Newton  # -initialThenCurrent

# Create the analysis object
analysis 	Static

# Run analysis
analyze 10