#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(graph
    PRIVATE
      DOF_Graph.cpp 
      Vertex.cpp 
      Graph.cpp
      DOF_GroupGraph.cpp  
      VertexIter.cpp
    PUBLIC
      DOF_Graph.h 
      Vertex.h 
      Graph.h
      DOF_GroupGraph.h  
      VertexIter.h
)


target_include_directories(graph PUBLIC ${CMAKE_CURRENT_LIST_DIR})
