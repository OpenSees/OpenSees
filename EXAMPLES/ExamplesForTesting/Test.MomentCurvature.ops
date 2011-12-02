# A procedure for performing section analysis (only does
# moment-curvature, but can be easily modified to do any mode
# of section reponse.
#
# MHS
# October 2000
#
# Arguments
#	secTag -- tag identifying section to be analyzed
#	axialLoad -- axial load applied to section (negative is compression)
#	maxK -- maximum curvature reached during analysis
#	numIncr -- number of increments used to reach maxK (default 100)
#
# Sets up a recorder which writes moment-curvature results to file
# section$secTag.out ... the moment is in column 1, and curvature in column 2

proc MomentCurvature {secTag axialLoad maxK {numIncr 100} } {
	# Define two nodes at (0,0)
	node 1 0.0 0.0
	node 2 0.0 0.0

	# Fix all degrees of freedom except axial and bending
	fix 1 1 1 1
	fix 2 0 1 0

	# Define element
	#                         tag ndI ndJ  secTag
	element zeroLengthSection  1   1   2  $secTag

	# Create recorder
	recorder Node section$secTag.out disp -time -node 2 -dof 3

	# Define constant axial load
	pattern Plain 1 "Constant" {
		load 2 $axialLoad 0.0 0.0
	}

	# Define analysis parameters
	integrator LoadControl 0 1 0 0
	system SparseGeneral -piv;	# Overkill, but may need the pivoting!
	test NormUnbalance 1.0e-9 10
	numberer Plain
	constraints Plain
	algorithm Newton
	analysis Static

	# Do one analysis for constant axial load
	analyze 1

	# Define reference moment
	pattern Plain 2 "Linear" {
		load 2 0.0 0.0 1.0
	}

	# Compute curvature increment
	set dK [expr $maxK/$numIncr]

	# Use displacement control at node 2 for section analysis
	integrator DisplacementControl 2 3 $dK 1 $dK $dK

	# Do the section analysis
	analyze $numIncr
}
