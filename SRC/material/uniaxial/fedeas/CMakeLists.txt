#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Uniaxial_Fedeas)
target_sources(OPS_Uniaxial_Fedeas
    PRIVATE
        ../FedeasMaterial.cpp
        FedeasBond1Material.cpp
        FedeasBond2Material.cpp
        FedeasConcr1Material.cpp
        FedeasConcr2Material.cpp
        FedeasConcr3Material.cpp
        FedeasHardeningMaterial.cpp
        FedeasHyster1Material.cpp
        FedeasHyster2Material.cpp
        FedeasSteel1Material.cpp
        FedeasSteel2Material.cpp
        PD1.c
        PlasticDamageMaterial.cpp
    PUBLIC
        ../FedeasMaterial.h
        FedeasBond1Material.h
        FedeasBond2Material.h
        FedeasConcr1Material.h
        FedeasConcr2Material.h
        FedeasConcr3Material.h
        FedeasHardeningMaterial.h
        FedeasHyster1Material.h
        FedeasHyster2Material.h
        FedeasSteel1Material.h
        FedeasSteel2Material.h
        PlasticDamageMaterial.h
)

add_library(OPS_Uniaxial_Fedeas_f OBJECT)
target_sources(OPS_Uniaxial_Fedeas_f 
    PRIVATE
        Hard1.f
        PD1.f
)
target_link_libraries(OPS_Uniaxial_Fedeas PUBLIC
  OPS_Uniaxial_Fedeas_f
  ${FEDEAS_Uniaxial_LIBRARIES}
  #${BUNDLE_LIBS}/fedeas.lib
)
#target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

