target_sources(libg3
    PRIVATE
      FEM_ObjectBroker.cpp
    PUBLIC
      FEM_ObjectBroker.h
)

if (OSP_OPTION_BROKER_ALL)
    target_sources(actor
      PRIVATE
        FEM_ObjectBrokerAllClasses.cpp
      PUBLIC
        FEM_ObjectBrokerAllClasses.h
)
endif()

target_include_directories(actor PUBLIC $(CMAKE_CURRENT_LIST_DIR))
