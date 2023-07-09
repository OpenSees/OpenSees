target_sources(OPS_Reliability
    PRIVATE
        NewmarkSensitivityIntegrator.cpp
        PFEMSensitivityIntegrator.cpp
        SensitivityAlgorithm.cpp
        SensitivityIntegrator.cpp
        StaticSensitivityIntegrator.cpp
    PUBLIC
        NewmarkSensitivityIntegrator.h
        PFEMSensitivityIntegrator.h
        SensitivityAlgorithm.h
        SensitivityIntegrator.h
        StaticSensitivityIntegrator.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

