#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_masonry)

target_sources(OPS_Element
  PUBLIC
  BeamGT.cpp
  MasonPan12.cpp
  MasonPan3D.cpp
  PRIVATE
  BeamGT.h
  MasonPan12.h
  MasonPan3D.h

)    

target_sources(OPS_Material
  PRIVATE
  "${OPS_SRC_DIR}/material/uniaxial/Masonry.cpp"
  "${OPS_SRC_DIR}/material/uniaxial/Masonryt.cpp"
  "${OPS_SRC_DIR}/material/uniaxial/Trilinwp.cpp"
  "${OPS_SRC_DIR}/material/uniaxial/Trilinwp2.cpp"
  "${OPS_SRC_DIR}/material/uniaxial/Trilinwpd.cpp"
  PUBLIC
  "${OPS_SRC_DIR}/material/uniaxial/Masonry.h"
  "${OPS_SRC_DIR}/material/uniaxial/Masonryt.h"
  "${OPS_SRC_DIR}/material/uniaxial/Trilinwp.h"
  "${OPS_SRC_DIR}/material/uniaxial/Trilinwp2.h"
  "${OPS_SRC_DIR}/material/uniaxial/Trilinwpd.h"
)
