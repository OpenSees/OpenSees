#!/home/fmk/bin/OpenSees

# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# Uniaxial Material Viewer
# ------------------------
#
# Written: fmk/MHS
# Date: March 2001

# ##############################################################
# Define some global variables & set their initial values
# ##############################################################

# current material id tag
set matID 0

# initial size of canvas used to draw stress-strain relations
set height 400
set width  400
set interpreterWidth 60
set interpreterHeight 8

# absolute max values for strain & stress drawn
set maxStress 50.0
set maxStrain  1.0

# last point plotted on canvas
set xLast [expr $height /2]
set yLast [expr  $width /2]

set minE -$maxStrain
set maxE $maxStrain

# some frames that will be toggled on & off
set toggleFrame .interpreter

# ##############################################################
# Set main window parameters
# ##############################################################

wm title . "Uniaxial Material"

# ##############################################################
# Invoke OpenSees command to make avaliable commands for testing
# the materials: 
#     1)uniaxialMaterial 
#     2)strainUniaxialTest 
#     3)stressUniaxialTest
#     4)tangUniaxialTest
# ##############################################################

model test 1

# ##############################################################
# Define menu frame at top of window 
# ##############################################################

frame .mbar -borderwidth 1 -relief raised
pack .mbar -fill x

menubutton .mbar.materials -text "Materials" -relief raised -menu .mbar.materials.menu
pack .mbar.materials -side left

button .mbar.values -text Values -command "SetValues"
pack .mbar.values -side left
button .mbar.settings -text Settings -command "Settings"
pack .mbar.settings -side left
button .mbar.reset -text Reset -command "Reset"
pack .mbar.reset -side left
button .mbar.quit -text Quit -command exit
pack .mbar.quit -side left

frame .style -borderwidth 1 -relief sunken
pack .style -fill x -pady 1

# ##############################################################
# Material Menu
# ##############################################################

set m [menu .mbar.materials.menu]

# source in the material data structures and procedures
source uniaxialMaterials.tcl

# ##############################################################
# Define the canvas and slider 
# ##############################################################

frame .figure 

# Define the scale

set theScale [ scale .figure.scale -from [expr -$maxStrain] -to $maxStrain \
	-length $width -variable strain -orient horizontal \
	-tickinterval [expr $maxStrain / 2] -resolution [expr $maxStrain / 100] \
	-showvalue true -command SetStrain ]

#pack .figure .figure.scale

# Define the canvas

set theCanvas [canvas .figure.canvas -bg #0088cc -height $height -width $width]

$theCanvas create line 10 [expr $height/2.0] [expr $width-10] [expr $height/2.0]
$theCanvas create line [expr $width/2.0] 10 [expr $width/2.0] [expr $height-10]
$theCanvas create text [expr $width-30] [expr $height/2.0-10] -text "Strain $maxStrain"
$theCanvas create text [expr $width/2.0+30] 10 -text "Stress $maxStress"

pack .figure .figure.canvas .figure.scale

pack .figure -side top

frame .style2 -borderwidth 1 -relief sunken
pack .style2 -fill x -pady 2

# ##############################################################
# Define the Settings frame & Settings procedure 
# ##############################################################

frame .settings -borderwidth 1 -relief raised

label .settings.lstrain -text "Strain (max absolute)" 
entry .settings.estrain -textvariable  maxStrain -relief sunken 
label .settings.lstress -text "Stress (max absolute)" 
entry .settings.estress -textvariable  maxStress -relief sunken 
button .settings.ok -text OK -command "Reset"

grid .settings.lstrain -row 0 -column 0 -sticky e
grid .settings.estrain -row 0 -column 1 -sticky ew
grid .settings.ok -row 0 -rowspan 2 -column 2 -sticky nsew
grid .settings.lstress -row 1 -column 0 -sticky e
grid .settings.estress -row 1 -column 1 -sticky ew

# Procedure Settings - used to set bottom frame to be the settings frame
proc Settings { } {
    global toggleFrame
    global .settings

    pack forget $toggleFrame
    set toggleFrame .settings
    pack $toggleFrame -side bottom  -fill x -pady 4
}


# ##############################################################
# Define the Values frame & SetValues procedure 
# ##############################################################

frame .values -borderwidth 1 -relief raised

label .values.strain -text "Strain : " 
label .values.stress -text "Stress : " 
label .values.tangent   -text "Tangent: " 
label .values.dstrain -text " 0.0" 
label .values.dstress -text " 0.0" 
label .values.dtangent   -text " 0.0" 

grid .values.strain -row 0 -column 0 -sticky e
grid .values.dstrain -row 0 -column 1 -sticky ew
grid .values.stress -row 1 -column 0 -sticky e
grid .values.dstress -row 1 -column 1 -sticky ew
grid .values.tangent -row 2 -column 0 -sticky e
grid .values.dtangent -row 2 -column 1 -sticky ew

# Procedure Settings - used to set bottom frame to be the settings frame
proc SetValues { } {
    global toggleFrame
    global .values

    pack forget $toggleFrame
    set toggleFrame .values
    pack $toggleFrame -side bottom  -fill x -pady 4
}

# ##############################################################
# Define the Reset Procedure
# ##############################################################

proc Reset { } {
    global theCanvas
    global theScale
    global height
    global width
    global maxStrain
    global maxStress
    global strain
    global xLast
    global yLast
    global width

    set strain 0

    $theScale configure -from [expr -$maxStrain] -to $maxStrain -length $width \
	    -tickinterval [expr $maxStrain / 2] -resolution [expr $maxStrain / 100] \
	    -showvalue true -command SetStrain 

    $theCanvas move all 400 400

    $theCanvas create line 10 [expr $height/2.0] [expr $width-10] [expr $height/2.0]
    $theCanvas create line [expr $width/2.0] 10 [expr $width/2.0] [expr $height-10]
    $theCanvas create text [expr $width-30] [expr $height/2.0-10] -text "Strain $maxStrain"
    $theCanvas create text [expr $width/2.0+30] 10 -text "Stress $maxStress"

    set xLast [expr $width/2]
    set yLast [expr $height/2]

    .values.dstrain  config -text [format "%6.4e" 0.0]
    .values.dstress  config -text [format "%6.4e" 0.0]
    .values.dtangent config -text [format "%6.4e" 0.0]

    ResetMaterial
}


proc ResetMaterial { } {
    global matID
    if {$matID != 0} {eval uniaxialTest $matID}
}


# ##############################################################
# Define the SetStrain Procedure
# ##############################################################

proc SetStrain {strain} {
    global theCanvas
    global maxStrain
    global maxStress
    global xLast
    global yLast
    global matID
    global height
    global width
    global toggleFrame
    
    if {$matID != 0} {
	eval strainUniaxialTest $strain
	set stress [stressUniaxialTest]
	
	set diffStrain [expr $width/(2*$maxStrain)]
	set diffStress [expr $height/(2*$maxStress)]

	set x [expr $width / 2 + $strain * $diffStrain]
	set y [expr $height / 2 - $stress * $diffStress]
	$theCanvas create line $xLast $yLast $x $y
	set xLast $x
	set yLast $y

	if {$toggleFrame == ".values"} {
	    set tangent [tangUniaxialTest]
	    .values.dstrain  config -text [format "%6.4e" $strain]
	    .values.dstress  config -text [format "%6.4e" $stress]
	    .values.dtangent config -text [format "%6.4e" $tangent]
	}
    }
}







