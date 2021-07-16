
target_sources(OPS_Material
    PRIVATE
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

add_library(OPS_Material_Uniaxial_Fedeas OBJECT)
target_sources(OPS_Material_Uniaxial_Fedeas 
    PRIVATE
        Hard1.f
        PD1.f
)
target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

