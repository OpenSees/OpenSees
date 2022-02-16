#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(matrix OBJECT)

target_sources(OPS_Matrix
    PRIVATE
      Matrix.cpp
      Vector.cpp
      ID.cpp
    PUBLIC
      Matrix.h
      Vector.h
      ID.h
)


target_include_directories(OPS_Matrix PUBLIC ${CMAKE_CURRENT_LIST_DIR})
