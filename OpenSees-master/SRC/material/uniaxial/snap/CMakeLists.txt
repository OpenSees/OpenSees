#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
    PRIVATE
    Pinching.cpp
    Bilinear.cpp
    Clough.cpp
    CloughHenry.cpp
    CloughDamage.cpp
    PinchingDamage.cpp
    PUBLIC
    Pinching.h
    Bilinear.h
    Clough.h
    CloughHenry.h
    CloughDamage.h
    PinchingDamage.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})


