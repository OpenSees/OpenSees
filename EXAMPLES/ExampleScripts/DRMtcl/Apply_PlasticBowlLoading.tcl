#================================================================
# Apply_LateralLoad.tcl
# --Create lateral load at pile head
# Zhaohui Yang & Boris Jeremic UC Davis
# Nov. 1, 2002
#================================================================

# Apply plastic bowl loading at specified elements using input accel. and displ. 
# and corresponding +x, -x,  +y, -y, +z and -z coordinates to specify the plastic bowl box
#Note: accel. and displ. are specified at all nodes
pattern PBowlLoading 1 -pbele "$Dir/PBElements.dat" -acce "$Dir/Inp_acce.dat" -disp "$Dir/Inp_disp.dat" -dt 0.02 -factor 1 -xp 6.0 -xm -6.0 -yp 6.0 -ym -6.0 -zp 0.0 -zm -17.5
					      
puts "Finished applying plastic bowl loading to specified elements..."
