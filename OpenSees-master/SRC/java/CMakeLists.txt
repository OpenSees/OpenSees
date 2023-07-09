#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Java OBJECT)

target_sources(OPS_Java
    PRIVATE
        OpenSeesEvaluator.cpp
    PUBLIC
        OpenSeesEvaluator.h
)
target_include_directories(OPS_Java PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(./kepler/)
