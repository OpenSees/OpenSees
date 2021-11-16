target_sources(OPS_Reliability
    PRIVATE
        DeltaFilter.cpp
        Filter.cpp
        FilterIter.cpp
        KooFilter.cpp
        StandardLinearOscillatorAccelerationFilter.cpp
        StandardLinearOscillatorDisplacementFilter.cpp
        StandardLinearOscillatorVelocityFilter.cpp
    PUBLIC
        DeltaFilter.h
        Filter.h
        FilterIter.h
        KooFilter.h
        StandardLinearOscillatorAccelerationFilter.h
        StandardLinearOscillatorDisplacementFilter.h
        StandardLinearOscillatorVelocityFilter.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

