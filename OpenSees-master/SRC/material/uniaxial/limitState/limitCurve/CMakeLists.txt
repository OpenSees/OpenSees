#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
    PRIVATE
        AxialCurve.cpp
        LimitCurve.cpp
        RotationShearCurve.cpp
        ShearCurve.cpp
        ThreePointCurve.cpp
        WrapperLimitCurve.cpp
    PUBLIC
        AxialCurve.h
        LimitCurve.h
        RotationShearCurve.h
        ShearCurve.h
        ThreePointCurve.h
        WrapperLimitCurve.h
)
target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

