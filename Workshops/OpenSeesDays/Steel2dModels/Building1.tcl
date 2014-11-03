# written: fmk
# units: kip & in

# set some variables
set gMotion el_centro
set scale 1.

set roofWeight  [expr 80*120.*72./1000.]; #kips; 
set floorWeight [expr 95*120.*72./1000.];

set percentLoadFrame [expr 15./120.] 
set percentMassFrame 0.5; # %mass frame takes

set dampRatio 0.03
set mode1 1
set mode2 3

set Fy 60.
set E 30000.
set b 0.03

# set up my lists
set floorOffsets {216. 150. 150. 150. 150. 150.}
set colOffsets   {288. 288. 288.} 
set colSizes     {W30X173 W30X173 W27X146 W27X146 W24X104 W24X104};
set colExtSizes  {W14X193 W14X193 W14X159 W14X159 W14X109 W14X109};
set beamSizes    {W30X99 W30X99 W27X94 W27X94 W24X76 W24X76};

source SteelMomentFrame2d_UniformExcitation.tcl


