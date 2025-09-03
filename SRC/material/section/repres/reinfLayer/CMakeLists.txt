#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Section_Repres
    PRIVATE
        CircReinfLayer.cpp
        RectReinfLayer.cpp	
        ReinfLayer.cpp
        StraightReinfLayer.cpp
    PUBLIC
        CircReinfLayer.h
        RectReinfLayer.h    
        ReinfLayer.h
        StraightReinfLayer.h
)
target_include_directories(OPS_Section_Repres PUBLIC ${CMAKE_CURRENT_LIST_DIR})

