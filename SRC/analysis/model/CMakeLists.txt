#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      AnalysisModel.cpp
      FE_EleIter.cpp
      DOF_GrpIter.cpp
    PUBLIC
      AnalysisModel.h
      FE_EleIter.h
      DOF_GrpIter.h
)

#target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})
