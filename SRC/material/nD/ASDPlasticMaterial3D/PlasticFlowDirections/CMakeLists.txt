#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PUBLIC
      AllPlasticFlowDirections.h
      VonMises_PF.h
      DruckerPrager_PF.h
      ConstantDilatancy_PF.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
