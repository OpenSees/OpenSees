#add_library(OPS_Element OBJECT)

target_sources(G3_Element
    PRIVATE
      Element.cpp
      ElementalLoad.cpp
      Information.cpp
    PUBLIC
      Element.h
      ElementalLoad.h
      Information.h
)

target_sources(OPS_Tcl
    PRIVATE
        TclElementCommands.cpp
    #PUBLIC
)

target_include_directories(OPS_Element PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(elasticBeamColumn)
add_subdirectory(zeroLength)

