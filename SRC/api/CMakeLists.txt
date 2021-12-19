#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Api OBJECT)

target_sources(OPS_InterpTcl
  PRIVATE
    elementAPI_TCL.cpp
)

target_sources(OPS_Api
    PUBLIC
        #elementAPI_Dummy.cpp
        packages.cpp
        #win32Functions.cpp
    PUBLIC
        elementAPI.h
        packages.h
)
target_include_directories(OPS_Api PUBLIC ${CMAKE_CURRENT_LIST_DIR})

