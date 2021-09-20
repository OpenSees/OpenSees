#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Section_Repres
    PRIVATE
        FiberSectionRepr.cpp
        SectionRepres.cpp
    PUBLIC
        FiberSectionRepr.h
        SectionRepres.h
)
target_include_directories(OPS_Section_Repres PUBLIC ${CMAKE_CURRENT_LIST_DIR})

