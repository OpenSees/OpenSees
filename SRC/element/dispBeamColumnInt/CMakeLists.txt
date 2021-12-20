#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_dispBeamColumnInt OBJECT)

target_sources(OPS_Element
    PRIVATE
        DispBeamColumn2dInt.cpp
        FiberSection2dInt.cpp
        LinearCrdTransf2dInt.cpp
        #TclDispBeamColumnIntCommand.cpp
    PUBLIC
        DispBeamColumn2dInt.h
        FiberSection2dInt.h
        LinearCrdTransf2dInt.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

