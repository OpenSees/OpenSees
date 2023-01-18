#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        ElastomericBearingBoucWen2d.cpp
        ElastomericBearingBoucWen3d.cpp
        ElastomericBearingBoucWenMod3d.cpp
        ElastomericBearingPlasticity2d.cpp
        ElastomericBearingPlasticity3d.cpp
        ElastomericBearingUFRP2d.cpp
        ElastomericX.cpp
        HDR.cpp
        LeadRubberX.cpp
        #TclElastomericBearingBoucWenCommand.cpp
        #TclElastomericBearingPlasticityCommand.cpp
        #TclElastomericBearingUFRPCommand.cpp
    PUBLIC
        ElastomericBearingBoucWen2d.h
        ElastomericBearingBoucWen3d.h
        ElastomericBearingBoucWenMod3d.h
        ElastomericBearingPlasticity2d.h
        ElastomericBearingPlasticity3d.h
        ElastomericBearingUFRP2d.h
        ElastomericX.h
        HDR.h
        LeadRubberX.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

