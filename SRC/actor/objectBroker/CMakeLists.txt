target_sources(G3_SysOfEqn
    PRIVATE
      FEM_ObjectBroker.cpp
    PUBLIC
      FEM_ObjectBroker.h
)

if (OSP_OPTION_BROKER_ALL)
    target_sources(G3_Actor
      PRIVATE
        FEM_ObjectBrokerAllClasses.cpp
      PUBLIC
        FEM_ObjectBrokerAllClasses.h
)
endif()

target_include_directories(G3_Actor PUBLIC $(CMAKE_CURRENT_LIST_DIR))
