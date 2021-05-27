
target_sources(graph
    PRIVATE
      RCM.cpp
      AMDNumberer.cpp
      SimpleNumberer.cpp
      GraphNumberer.cpp
    PUBLIC
      RCM.h
      AMDNumberer.h
      SimpleNumberer.h
      GraphNumberer.h
)


target_include_directories(graph PUBLIC $(CMAKE_CURRENT_LIST_DIR))
