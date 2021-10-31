#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(G3_ObjectBroker
    PRIVATE
      FEM_ObjectBroker.cpp
    PUBLIC
      FEM_ObjectBroker.h
)

target_sources(OPS_ObjectBroker
    PRIVATE
      FEM_ObjectBrokerAllClasses.cpp
    PUBLIC
      FEM_ObjectBrokerAllClasses.h
)
target_link_libraries(OPS_ObjectBroker PUBLIC G3_ObjectBroker)
target_include_directories(OPS_ObjectBroker PUBLIC ${CMAKE_CURRENT_LIST_DIR})
