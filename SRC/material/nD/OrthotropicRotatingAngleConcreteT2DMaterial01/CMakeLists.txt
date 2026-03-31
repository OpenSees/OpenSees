#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
      OrthotropicRotatingAngleConcreteT2DMaterial01.cpp
    PUBLIC
      OrthotropicRotatingAngleConcreteT2DMaterial01.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
