#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_InterpTcl
    PRIVATE
        myCommands.cpp
)

target_sources(OPS_INTERPRETER
    PRIVATE
        Block2D.cpp
        Block3D.cpp
    PUBLIC
        Block2D.h
        Block3D.h
)

target_include_directories(OPS_InterpTcl PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_INTERPRETER PUBLIC ${CMAKE_CURRENT_LIST_DIR})

