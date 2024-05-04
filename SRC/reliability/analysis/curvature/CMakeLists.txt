
target_sources(OPS_Reliability
    PRIVATE
        CurvatureFitting.cpp
        CurvaturesBySearchAlgorithm.cpp
        FindCurvatures.cpp
        FirstPrincipalCurvature.cpp
    PUBLIC
        CurvatureFitting.h
        CurvaturesBySearchAlgorithm.h
        FindCurvatures.h
        FirstPrincipalCurvature.h
)
target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})

