target_sources(OPS_Actor
    PRIVATE
      FEM_ObjectBroker.cpp
      FEM_ObjectBrokerAllClasses.cpp
    PUBLIC
      FEM_ObjectBroker.h
      FEM_ObjectBrokerAllClasses.h
)

target_include_directories(OPS_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})
