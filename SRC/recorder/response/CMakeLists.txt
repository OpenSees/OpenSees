#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Recorder
    PRIVATE 
        Response.cpp
        MaterialResponse.cpp
        ElementResponse.cpp
        FiberResponse.cpp
        CompositeResponse.cpp
    PUBLIC 
        Response.h
        MaterialResponse.h
        ElementResponse.h
        FiberResponse.h
        CompositeResponse.h
)

