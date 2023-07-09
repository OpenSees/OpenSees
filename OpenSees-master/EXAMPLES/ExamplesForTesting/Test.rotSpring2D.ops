# rotSpring2D.tcl
# Procedure which creates a rotational spring for a planar problem
#
# SETS A MULTIPOINT CONSTRAINT ON THE TRANSLATIONAL DEGREES OF FREEDOM,
# SO DO NOT USE THIS PROCEDURE IF THERE ARE TRANSLATIONAL ZEROLENGTH
# ELEMENTS ALSO BEING USED BETWEEN THESE TWO NODES
#
# Written: MHS
# Date: Jan 2000
#
# Formal arguments
#	eleID - unique element ID for this zero length rotational spring
#	nodeR - node ID which will be retained by the multi-point constraint
#	nodeC - node ID which will be constrained by the multi-point constraint
#	matID - material ID which represents the moment-rotation relationship
#		for the spring

proc rotSpring2D {eleID nodeR nodeC matID} {
	# Create the zero length element
	element zeroLength $eleID $nodeR $nodeC -mat $matID -dir 6

	# Constrain the translational DOF with a multi-point constraint
	#          retained constrained DOF_1 DOF_2 ... DOF_n
	equalDOF    $nodeR     $nodeC     1     2
}