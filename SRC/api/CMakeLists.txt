add_library(OPS_Api OBJECT)

target_sources(OPS_Api
    PRIVATE
        #elementAPI_Dummy.cpp
        elementAPI_TCL.cpp
        packages.cpp
	    #win32Functions.cpp
    PUBLIC
        elementAPI.h
        packages.h
)
target_include_directories(OPS_Api PUBLIC ${CMAKE_CURRENT_LIST_DIR})

