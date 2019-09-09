proc DisplayPlane {ShapeType dAmp viewPlane {nEigen 0}  {quadrant 0}} {
	######################################################################################
	## DisplayPlane $ShapeType $dAmp $viewPlane $nEigen $quadrant
	######################################################################################
	## setup display parameters for specified viewPlane and display
	## 			Silvia Mazzoni & Frank McKenna, 2006
	##
	## 	ShapeType : 	type of shape to display. # options: ModeShape , NodeNumbers , DeformedShape 
	## 	dAmp : 		relative amplification factor for deformations
	## 	viewPlane :	set local xy axes in global coordinates (XY,YX,XZ,ZX,YZ,ZY)
	## 	nEigen : 		if nEigen not=0, show mode shape for nEigen eigenvalue
	##	quadrant:		quadrant where to show this figure (0=full figure)
	##	
	######################################################################################

	set Xmin [lindex [nodeBounds] 0];	# view bounds in global coords -  will add padding on the sides
	set Ymin [lindex [nodeBounds] 1];
	set Zmin [lindex [nodeBounds] 2];
	set Xmax [lindex [nodeBounds] 3];
	set Ymax [lindex [nodeBounds] 4];
	set Zmax [lindex [nodeBounds] 5];

	set Xo 0;	# center of local viewing system
	set Yo 0;
	set Zo 0;

	set uLocal [string index $viewPlane 0];	# viewPlane local-x axis in global coordinates
	set vLocal [string index $viewPlane 1];	# viewPlane local-y axis in global coordinates


	if  {$viewPlane =="3D" } {
		set uMin $Zmin+$Xmin
		set uMax $Zmax+$Xmax
		set vMin $Ymin
		set vMax $Ymax
		set wMin -10000
		set wMax 10000
		vup 0 1 0; # dirn defining up direction of view plane
	} else {
		set keyAxisMin "X $Xmin Y $Ymin Z $Zmin"
		set keyAxisMax "X $Xmax Y $Ymax Z $Zmax"
		set axisU [string index $viewPlane 0];
		set axisV [string index $viewPlane 1];
		set uMin [string map $keyAxisMin $axisU]
		set uMax [string map $keyAxisMax $axisU]
		set vMin [string map $keyAxisMin $axisV]
		set vMax [string map $keyAxisMax $axisV]
		if {$viewPlane =="YZ" || $viewPlane =="ZY" } {
			set wMin $Xmin
			set wMax $Xmax
		} elseif  {$viewPlane =="XY" || $viewPlane =="YX" } {
			set wMin $Zmin
			set wMax $Zmax
		} elseif  {$viewPlane =="XZ" || $viewPlane =="ZX" } {
			set wMin $Ymin
			set wMax $Ymax
		} else {
		return -1
		}
	}

	set epsilon 1e-3;	# make windows width or height not zero when the Max and Min values of a coordinate are the same

	set uWide [expr $uMax - $uMin+$epsilon];
	set vWide [expr $vMax - $vMin+$epsilon];
	set uSide [expr 0.25*$uWide];
	set vSide [expr 0.25*$vWide];
	set uMin [expr $uMin - $uSide];
	set uMax [expr $uMax + $uSide];
	set vMin [expr $vMin - $vSide];
	set vMax [expr $vMax + 2*$vSide];	# pad a little more on top, because of window title
	set uWide [expr $uMax - $uMin+$epsilon];
	set vWide [expr $vMax - $vMin+$epsilon];
	set uMid [expr ($uMin+$uMax)/2];
	set vMid [expr ($vMin+$vMax)/2];

	# keep the following general, as change the X and Y and Z for each view plane
	# next three commmands define viewing system, all values in global coords
	vrp $Xo $Yo $Zo;    # point on the view plane in global coord, center of local viewing system
	if {$vLocal == "X"} {
		vup 1 0 0; # dirn defining up direction of view plane
	} elseif {$vLocal == "Y"} {
		vup 0 1 0; # dirn defining up direction of view plane
	} elseif {$vLocal == "Z"} {
		vup 0 0 1; # dirn defining up direction of view plane
	}
	if {$viewPlane =="YZ" } {
		vpn 1 0 0; # direction of outward normal to view plane
		prp 10000. $uMid $vMid ; # eye location in local coord sys defined by viewing system
		plane 10000 -10000; # distance to front and back clipping planes from eye
	} elseif  {$viewPlane =="ZY" } {
		vpn -1 0 0; # direction of outward normal to view plane
		prp -10000. $vMid $uMid ; # eye location in local coord sys defined by viewing system
		plane 10000 -10000; # distance to front and back clipping planes from eye
	} elseif  {$viewPlane =="XY"  } {
		vpn 0 0 1; # direction of outward normal to view plane
		prp $uMid $vMid 10000; # eye location in local coord sys defined by viewing system
		plane 10000 -10000; # distance to front and back clipping planes from eye
	} elseif  {$viewPlane =="YX" } {
		vpn 0 0 -1; # direction of outward normal to view plane
		prp $uMid $vMid -10000; # eye location in local coord sys defined by viewing system
		plane 10000 -10000; # distance to front and back clipping planes from eye
	} elseif  {$viewPlane =="XZ" } {
		vpn 0 -1 0; # direction of outward normal to view plane
		prp $uMid -10000 $vMid ; # eye location in local coord sys defined by viewing system
		plane 10000 -10000; # distance to front and back clipping planes from eye
	} elseif  {$viewPlane =="ZX" } {
		vpn 0 1 0; # direction of outward normal to view plane
		prp $uMid 10000 $vMid ; # eye location in local coord sys defined by viewing system
		plane 10000 -10000; # distance to front and back clipping planes from eye
	} elseif  {$viewPlane =="3D" } {
		vpn 1 0.25 1.25; # direction of outward normal to view plane
		prp -100 $vMid 10000; # eye location in local coord sys defined by viewing system
		plane 10000 -10000; # distance to front and back clipping planes from eye
	}  else {
		return -1
	}
	# next three commands define view, all values in local coord system
	if  {$viewPlane =="3D" } {
		viewWindow [expr $uMin-$uWide/4] [expr $uMax/2] [expr $vMin-0.25*$vWide] [expr $vMax] 
	} else {
		viewWindow $uMin $uMax $vMin $vMax
	}
	projection 1; 	# projection mode, 0:prespective, 1: parallel
	fill 1; 		# fill mode; needed only for solid elements

	if {$quadrant == 0} {
		port -1 1 -1 1 	# area of window that will be drawn into (uMin,uMax,vMin,vMax);
	} elseif {$quadrant == 1} {
		port 0 1 0 1 	# area of window that will be drawn into (uMin,uMax,vMin,vMax);
	} elseif {$quadrant == 2} {
		port -1 0 0 1 	# area of window that will be drawn into (uMin,uMax,vMin,vMax);
	} elseif {$quadrant == 3} {
		port -1 0 -1 0 	# area of window that will be drawn into (uMin,uMax,vMin,vMax);
	} elseif {$quadrant == 4} {
		port 0 1 -1 0 	# area of window that will be drawn into (uMin,uMax,vMin,vMax);
	}

	if {$ShapeType ==  "ModeShape" } {
		display -$nEigen 0  [expr 5.*$dAmp]; 	# display mode shape for mode $nEigen
	} elseif  {$ShapeType ==  "NodeNumbers" } {
		display 1 -1 0  ; 		# display node numbers
	} elseif  {$ShapeType ==  "DeformedShape" }  {
		display 1 2 $dAmp; 		# display deformed shape  the 2 makes the nodes small
	}
};                                                                                                                                                          #
######################################################################################

