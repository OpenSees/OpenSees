
target_sources(OPS_Reliability
    PRIVATE
        Hessian.cpp
        HessianApproximation.cpp
    PUBLIC
        Hessian.h
        HessianApproximation.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

