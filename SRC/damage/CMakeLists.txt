#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Damage
  PRIVATE
    DamageModel.cpp
    HystereticEnergy.cpp
    Kratzig.cpp
    Mehanny.cpp
    ParkAng.cpp
    NormalizedPeak.cpp
    DamageResponse.cpp
  PUBLIC
    DamageModel.h
    HystereticEnergy.h
    Kratzig.h
    Mehanny.h
    ParkAng.h
    NormalizedPeak.h
    DamageResponse.h
)

target_include_directories(OPS_Damage PUBLIC ${CMAKE_CURRENT_LIST_DIR})

