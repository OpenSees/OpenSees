#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      Analysis.cpp 
      DirectIntegrationAnalysis.cpp 
      DomainDecompositionAnalysis.cpp
      DomainUser.cpp 
      EigenAnalysis.cpp
      ResponseSpectrumAnalysis.cpp
      SDFAnalysis.cpp
      StaticAnalysis.cpp 
      StaticDomainDecompositionAnalysis.cpp 
      SubstructuringAnalysis.cpp    
      TransientAnalysis.cpp
      TransientDomainDecompositionAnalysis.cpp 
      VariableTimeStepDirectIntegrationAnalysis.cpp
    PUBLIC
      Analysis.h 
      DirectIntegrationAnalysis.h 
      DomainDecompositionAnalysis.h
      DomainUser.h 
      EigenAnalysis.h
      ResponseSpectrumAnalysis.h
      StaticAnalysis.h 
      StaticDomainDecompositionAnalysis.h 
      SubstructuringAnalysis.h    
      TransientAnalysis.h
      TransientDomainDecompositionAnalysis.h 
      VariableTimeStepDirectIntegrationAnalysis.h
)
target_sources(OPS_PFEM
    PRIVATE
        PFEMAnalysis.cpp 
    PUBLIC
        PFEMAnalysis.h 
)
#target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#add_subdirectory(analysis)
