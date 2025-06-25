#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_library(EIGENAPI INTERFACE)

target_include_directories(EIGENAPI INTERFACE ${CMAKE_CURRENT_LIST_DIR})

target_sources(EIGENAPI
    INTERFACE
      EigenAPI.h
      converters.h
      operations.h
      typedefs.h
)
