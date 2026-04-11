#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Section_Repres
    PRIVATE
        CircPatch.cpp
        QuadPatch.cpp
    PUBLIC
        CircPatch.h
        Patch.h
        QuadPatch.h
)
target_include_directories(OPS_Section_Repres PUBLIC ${CMAKE_CURRENT_LIST_DIR})

