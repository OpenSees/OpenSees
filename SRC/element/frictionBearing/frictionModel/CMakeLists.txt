
target_sources(OPS_Element_frictionBearing
    PRIVATE
        Coulomb.cpp
        FrictionModel.cpp
        FrictionResponse.cpp
        VelDependent.cpp
        VelDepMultiLinear.cpp
        VelNormalFrcDep.cpp
        VelPressureDep.cpp
    PUBLIC
        Coulomb.h
        FrictionModel.h
        FrictionResponse.h
        VelDependent.h
        VelDepMultiLinear.h
        VelNormalFrcDep.h
        VelPressureDep.h
)
target_include_directories(OPS_Element_frictionBearing PUBLIC ${CMAKE_CURRENT_LIST_DIR})

