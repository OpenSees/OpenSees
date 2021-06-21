target_sources(LG3_Analysis
    PRIVATE
      AnalysisModel.cpp
      FE_EleIter.cpp
      DOF_GrpIter.cpp
    PUBLIC
      AnalysisModel.h
      FE_EleIter.h
      DOF_GrpIter.h
)

#target_include_directories(LG3_Analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
