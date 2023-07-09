#----------------------------------------------------------
#  create elastic pile section
#----------------------------------------------------------

section Elastic 1  25000000 0.785 0.049 0.049 9615385 0.098

# elastic torsional material for combined 3D section
uniaxialMaterial Elastic 3000   1.e10

# create combined 3D section
set secTag3D 3
section Aggregator $secTag3D  3000   T   -section 1
