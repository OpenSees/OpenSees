#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
  PRIVATE
    DrainHardeningMaterial.cpp
    DrainBilinearMaterial.cpp
    DrainClough1Material.cpp
    DrainClough2Material.cpp
    DrainPinch1Material.cpp
  PUBLIC
    DrainHardeningMaterial.h
    DrainBilinearMaterial.h
    DrainClough1Material.h
    DrainClough2Material.h
    DrainPinch1Material.h
)

add_library(OPS_Material_Uniaxial_Drain_f)
target_sources(OPS_Material_Uniaxial_Drain_f PRIVATE anal00.f common00.f)
target_link_libraries(OPS_Material PUBLIC OPS_Material_Uniaxial_Drain_f) 

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})


