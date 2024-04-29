#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
      # ASDPlasticMaterialGlobals.cpp
      OPS_AllASDPlasticMaterials.cpp
    PUBLIC
      AllASDPlasticMaterials.h     
      ASDPlasticMaterial.h        
      ElasticityBase.h    
      YieldFunctionBase.h
      ASDPlasticMaterialGlobals.h  
      ASDPlasticMaterialTraits.h  
      PlasticFlowBase.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(ElasticityModels)
add_subdirectory(PlasticFlowDirections)
add_subdirectory(YieldFunctions)
