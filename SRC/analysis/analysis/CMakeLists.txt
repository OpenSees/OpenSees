
target_sources(analysis
    PRIVATE
      DomainUser.cpp 
      Analysis.cpp 
      StaticAnalysis.cpp 
      TransientAnalysis.cpp
      DirectIntegrationAnalysis.cpp 
      DomainDecompositionAnalysis.cpp
      SubstructuringAnalysis.cpp	
      EigenAnalysis.cpp
      VariableTimeStepDirectIntegrationAnalysis.cpp
      StaticDomainDecompositionAnalysis.cpp 
      TransientDomainDecompositionAnalysis.cpp 
      PFEMAnalysis.cpp 
    PUBLIC
      DomainUser.h 
      Analysis.h 
      StaticAnalysis.h 
      TransientAnalysis.h
      DirectIntegrationAnalysis.h 
      DomainDecompositionAnalysis.h
      SubstructuringAnalysis.h	
      EigenAnalysis.h
      VariableTimeStepDirectIntegrationAnalysis.h
      StaticDomainDecompositionAnalysis.h 
      TransientDomainDecompositionAnalysis.h 
      PFEMAnalysis.h 
)

target_include_directories(analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))

#add_subdirectory(analysis)
