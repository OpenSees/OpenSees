# a window showing the displaced shape
#recorder display g3 10 10 300 300 -wipe

# next three commands define viewing system, all values in global coords
vrp 288.0 150.0 0    # point on the view plane in global coord, center of local viewing system
vup 0 1 0            # dirn defining up direction of view plane
vpn 0 0 1            # direction of outward normal to view plane

# next three commands define view, all values in local coord system
prp 0 0 100                   # eye location in local coord sys defined by viewing system
viewWindow -400 400 -400 400  # view bounds uMin, uMax, vMin, vMax in local coords
plane 0 150                   # distance to front and back clipping planes from eye
projection 0                  # projection mode

port -1 1 -1 1                # area of window that will be drawn into
fill 1                        # fill mode
display 1 0 10                


