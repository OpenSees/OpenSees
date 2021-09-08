target_sources(actor
    PRIVATE
      FEM_ObjectBroker.cpp
      FEM_ObjectBrokerAllClasses.cpp
    PUBLIC
      FEM_ObjectBroker.h
      FEM_ObjectBrokerAllClasses.h
)

target_include_directories(actor PUBLIC $(CMAKE_CURRENT_LIST_DIR))
