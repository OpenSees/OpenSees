#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material_YieldSurface
    PRIVATE
        ExponReducing.cpp
        MultiLinearKp.cpp
        NullPlasticMaterial.cpp
        PlasticHardeningMaterial.cpp
    PUBLIC
        ExponReducing.h
        MultiLinearKp.h
        NullPlasticMaterial.h
        PlasticHardeningMaterial.h
)
target_include_directories(OPS_Material_YieldSurface PUBLIC ${CMAKE_CURRENT_LIST_DIR})

