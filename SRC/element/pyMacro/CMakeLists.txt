#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_pyMacro OBJECT)

target_sources(OPS_Element
    PRIVATE
        PY_Macro2D.cpp
    PUBLIC
        PY_Macro2D.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

