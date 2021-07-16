add_library(OPS_Element_generic OBJECT)

target_sources(OPS_Element_generic
    PRIVATE
        GenericClient.cpp
        GenericCopy.cpp
        #TclGenericClientCommand.cpp
        #TclGenericCopyCommand.cpp
    PUBLIC
        GenericClient.h
        GenericCopy.h
)
target_include_directories(OPS_Element_generic PUBLIC ${CMAKE_CURRENT_LIST_DIR})

