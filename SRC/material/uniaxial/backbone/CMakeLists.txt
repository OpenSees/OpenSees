#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
    PRIVATE
        ArctangentBackbone.cpp
        CappedBackbone.cpp
        HystereticBackbone.cpp
        LinearCappedBackbone.cpp
        ManderBackbone.cpp
        MaterialBackbone.cpp
        MultilinearBackbone.cpp
        RaynorBackbone.cpp
        ReeseSandBackbone.cpp
        ReeseSoftClayBackbone.cpp
        ReeseStiffClayBelowWS.cpp
        TrilinearBackbone.cpp
    PUBLIC
        ArctangentBackbone.h
        CappedBackbone.h
        HystereticBackbone.h
        LinearCappedBackbone.h
        ManderBackbone.h
        MaterialBackbone.h
        MultilinearBackbone.h
        RaynorBackbone.h
        ReeseSandBackbone.h
        ReeseSoftClayBackbone.h
        ReeseStiffClayBelowWS.h
        TrilinearBackbone.h
)
target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

