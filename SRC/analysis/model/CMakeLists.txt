target_sources(analysis
    PRIVATE
      AnalysisModel.cpp
      FE_EleIter.cpp
      DOF_GrpIter.cpp
    PUBLIC
      AnalysisModel.h
      FE_EleIter.h
      DOF_GrpIter.h
)

target_include_directories(analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
