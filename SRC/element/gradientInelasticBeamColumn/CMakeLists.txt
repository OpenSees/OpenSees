#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

#add_library(OPS_Element_gradientInelasticBeamColumn OBJECT)

target_sources(OPS_Element
    PRIVATE
        GradientInelasticBeamColumn2d.cpp
        GradientInelasticBeamColumn3d.cpp
        #TclGradientInelasticBeamColumnCommand.cpp
    PUBLIC
        GradientInelasticBeamColumn2d.h
        GradientInelasticBeamColumn3d.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

