add_library(element "")

target_sources(element
    PRIVATE
      Element.cpp
      ElementalLoad.cpp
      Information.cpp
    PUBLIC
      Element.h
      ElementalLoad.h
      Information.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(elasticBeamColumn)
