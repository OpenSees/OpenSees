#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        MixedBeamColumn2d.cpp
        MixedBeamColumn3d.cpp
        MixedBeamColumnAsym3d.cpp
    PUBLIC
        MixedBeamColumn2d.h
        MixedBeamColumn3d.h
        MixedBeamColumnAsym3d.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

