######################################################################
# Main file for py pushover analysis to illustrate PySimple1Gen
# command
#
# Created by Scott Brandenberg, April 30, 2004.
#######################################################################

#############################################################
# BUILD MODEL
#
#create the ModelBuilder for the zeroLength elements
model basic -ndm 2 -ndf 3

# Use PySimple1Gen and TzSimple1Gen commands
PySimple1Gen "PySoilProp.txt" "Nodes.tcl" "PySimple1.tcl" "PileElements.tcl" "PyMaterials.tcl" "PyPattern.tcl"
TzSimple1Gen "TzSoilProp.txt" "Nodes.tcl" "TzSimple1.tcl" "PileElements.tcl" "TzMaterials.tcl"
