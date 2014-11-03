# written: fmk
# units: kip & in

# set some variables
set gMotion el_centro
set scale 1.

set roofWeight  0.; #kips; 
set floorWeight 0.;
 
set percentMassFrame  1; # %mass taken by frame for transient analysis
set percentLoadFrame  1; # %total load taken by frame for gravity analysis

set dampRatio 0.0
set mode1 1
set mode2 2

set Fy 60.
set E 30000.
set b 0.03

# set up my lists
set floorOffsets {}
set colOffsets   {} 
set colSizes     {};
set colExtSizes  {};
set beamSizes    {};

source SteelMomentFrame2d_UniformExcitation.tcl

