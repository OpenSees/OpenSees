#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
      SAniSandMS.cpp
      SAniSandMS3D.cpp
      SAniSandMSPlaneStrain.cpp
    PUBLIC
      SAniSandMS.h
      SAniSandMS3D.h
      SAniSandMSPlaneStrain.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
