#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_InterpTcl
  PRIVATE
    elementAPI_TCL.cpp
    packages.cpp
)

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")

add_library(OPS_Api OBJECT)

target_sources(OPS_Api
    PUBLIC
        #elementAPI_Dummy.cpp
        win32Functions.cpp
    PUBLIC
        elementAPI.h
        GenericDict.h
	#packages.h
)

target_include_directories(OPS_Api PUBLIC ${CMAKE_CURRENT_LIST_DIR})

endif()


