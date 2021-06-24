
target_sources(G3_SysOfEqn
    PRIVATE
        UmfpackGenLinSOE.cpp
        UmfpackGenLinSolver.cpp

    PUBLIC
        UmfpackGenLinSOE.h
        UmfpackGenLinSolver.h

)

target_include_directories(OPS_SysOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))

