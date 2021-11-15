#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_UP_ucsd OBJECT)

#target_sources(OPS_Element_UP_ucsd
target_sources(OPS_Element
    PRIVATE
        BBarBrickUP.cpp
        BBarFourNodeQuadUP.cpp
        BrickUP.cpp
        FourNodeQuadUP.cpp
        Nine_Four_Node_QuadUP.cpp
    #Nine_Four_Node_QuadUPOld.cpp
        shp3dv.cpp
        #TclFourNodeQuadUPCommand.cpp
        Twenty_Eight_Node_BrickUP.cpp
    PUBLIC
        BBarBrickUP.h
        BBarFourNodeQuadUP.h
        BrickUP.h
        FourNodeQuadUP.h
        Nine_Four_Node_QuadUP.h
    #Nine_Four_Node_QuadUPOld.h
        shp3dv.h
        Twenty_Eight_Node_BrickUP.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

