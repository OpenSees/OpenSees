#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_frictionBearing OBJECT)

target_sources(OPS_Element
    PRIVATE
        FlatSliderSimple2d.cpp
        FlatSliderSimple3d.cpp
        FPBearingPTV.cpp
        MultiFP2d.cpp
        RJWatsonEQS2d.cpp
        RJWatsonEQS3d.cpp
        SingleFPSimple2d.cpp
        SingleFPSimple3d.cpp
        #TclFlatSliderCommand.cpp
        #TclRJWatsonEQSCommand.cpp
        #TclSingleFPCommand.cpp
        TFP_Bearing2d.cpp
        TFP_Bearing.cpp
        TPB1D.cpp
        TripleFrictionPendulum.cpp
        TripleFrictionPendulumX.cpp	
    PUBLIC
        FlatSliderSimple2d.h
        FlatSliderSimple3d.h
        FPBearingPTV.h
        MultiFP2d.h
        RJWatsonEQS2d.h
        RJWatsonEQS3d.h
        SingleFPSimple2d.h
        SingleFPSimple3d.h
        TFP_Bearing2d.h
        TFP_Bearing.h
        TPB1D.h
        TripleFrictionPendulum.h
        TripleFrictionPendulumX.h	
)

add_subdirectory(frictionModel)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

