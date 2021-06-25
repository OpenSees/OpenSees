
target_sources(OPS_Tcl
    PRIVATE
        Block2D.cpp
        Block3D.cpp
        myCommands.cpp
    PUBLIC
        Block2D.h
        Block3D.h
)
target_include_directories(OPS_Tcl PUBLIC $(CMAKE_CURRENT_LIST_DIR))

