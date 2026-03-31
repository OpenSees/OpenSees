#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
      # ASDPlasticMaterial3DGlobals.cpp
      OPS_AllASDPlasticMaterial3Ds.cpp
    PUBLIC
      AllASDPlasticMaterial3Ds.h     
      ASDPlasticMaterial3D.h        
      ElasticityBase.h    
      YieldFunctionBase.h
      ASDPlasticMaterial3DGlobals.h  
      ASDPlasticMaterial3DTraits.h  
      PlasticFlowBase.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(ElasticityModels)
add_subdirectory(PlasticFlowDirections)
add_subdirectory(YieldFunctions)
