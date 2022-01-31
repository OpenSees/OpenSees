#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Element
    PRIVATE
    Brick.cpp
    Twenty_Node_Brick.cpp
    shp3d.cpp
    BbarBrick.cpp
    BbarBrickWithSensitivity.cpp
    PUBLIC
    Brick.h
    Twenty_Node_Brick.h
    shp3d.h
    BbarBrick.h
    BbarBrickWithSensitivity.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
