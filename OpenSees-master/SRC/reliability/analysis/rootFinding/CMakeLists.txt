
target_sources(OPS_Reliability
    PRIVATE
        # ModNewtonRootFinding.cpp
        RootFinding.cpp
        SafeGuardedZeroFindingAlgorithm.cpp
        SecantRootFinding.cpp
        ZeroFindingAlgorithm.cpp
    PUBLIC
        # ModNewtonRootFinding.h
        RootFinding.h
        SafeGuardedZeroFindingAlgorithm.h
        SecantRootFinding.h
        ZeroFindingAlgorithm.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

