# ************************************************************************************************************
#              ***************************************************************************
#                                  *************************************
#									Collaboration: Aryan, Petras, Iman                 
#                                  *************************************
#              ***************************************************************************
# ************************************************************************************************************


wipe

model	basic	-ndm	3	-ndf	6;

set			AnalysisStartT [clock seconds] 

source		01_Output/01_Output_Folder.tcl
puts		"-----------------------------"
puts		""
# ----------------------------------------------------------------
source		02_Parameters/02_00_Main_Paramemetrs.tcl
puts		"Main Parameters: Defied"
puts		"-----------------------------"
puts		""
# ----------------------------------------------------------------

source		02_Parameters/02_00_1_Plate_Parameters.tcl
puts		"Plate Parameters: Defied"
puts		"-----------------------------"
puts		""
# ----------------------------------------------------------------

source		04_Nodes/04_01_Nodes.tcl
puts		"Nodes: Defied"
puts		"-----------------------------"
puts		""
# ----------------------------------------------------------------

source		05_Elements/05_01_Transformation.tcl
source 		05_Elements/05_02_Elements.tcl
puts		"Elements: Defied"
puts		"-----------------------------"
puts		""
# ----------------------------------------------------------------

source		06_Node_Constraints/06_01.tcl
puts		"Boundary Conditions: Defied"
puts		"-----------------------------"
puts		""
# ----------------------------------------------------------------

source		07_Loads/07_01.tcl
puts		"Loads: Defied"
puts		"-----------------------------"
puts		""
# ----------------------------------------------------------------

source		08_Display_Modes/08_1_DisplayModel3D.tcl
source		08_Display_Modes/08_2_DisplayPlane.tcl
puts		"Display Modes: Defied"
puts		"-----------------------------"
puts		""
# ----------------------------------------------------------------

source		09_Recorders/09_Recorders.tcl
puts		"Recorders: Defied"
puts		"-----------------------------"
puts		""
# ----------------------------------------------------------------

puts		"Analysis: Started"
source		10_Analysis/10_Analysis.tcl
puts		"-----------------------------"
puts		""
# ----------------------------------------------------------------

puts 		"Full Analysis Done, Successfully"

set			AnalysisEndT 				[clock seconds];
set			AnalysisRunTime_Minute 		[expr 1.*($AnalysisEndT - $AnalysisStartT)/60.0];
puts		"Run Time = $AnalysisRunTime_Minute Minute(s)"
