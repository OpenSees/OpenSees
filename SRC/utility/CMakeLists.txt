#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Utilities
    PRIVATE
    Timer.cpp 
    FileIter.cpp 
    File.cpp 
    SimulationInformation.cpp 
    StringContainer.cpp
    PeerNGA.cpp
    PUBLIC
    Timer.h 
    FileIter.h 
    File.h 
    SimulationInformation.h 
    StringContainer.h 
)

target_include_directories(OPS_Utilities PUBLIC ${CMAKE_CURRENT_LIST_DIR})
