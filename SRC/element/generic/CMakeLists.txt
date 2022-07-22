#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        GenericClient.cpp
        GenericCopy.cpp
        #TclGenericClientCommand.cpp
        #TclGenericCopyCommand.cpp
    PUBLIC
        GenericClient.h
        GenericCopy.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

