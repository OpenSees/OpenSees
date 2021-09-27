#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
    PRIVATE
        LimitStateMaterial.cpp
        PinchingLimitStateMaterial.cpp
    PUBLIC
        LimitStateMaterial.h
        PinchingLimitStateMaterial.h
)
target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(./limitCurve/)

