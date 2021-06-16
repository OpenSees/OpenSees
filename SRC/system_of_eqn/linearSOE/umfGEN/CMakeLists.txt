
target_sources(systemOfEqn
    PRIVATE
        UmfpackGenLinSOE.cpp
        UmfpackGenLinSolver.cpp

    PUBLIC
        UmfpackGenLinSOE.h
        UmfpackGenLinSolver.h

)

target_include_directories(systemOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))

