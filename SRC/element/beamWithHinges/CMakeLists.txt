#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_library(OPS_Element_beamWithHinges OBJECT)

target_sources(OPS_Element_beamWithHinges
    PRIVATE
        #TclBeamWithHingesBuilder.cpp
    PUBLIC
)

target_include_directories(OPS_Element_beamWithHinges PUBLIC ${CMAKE_CURRENT_LIST_DIR})

